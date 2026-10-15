#pragma once
#include "../IterativeSolver.h"
#include "BlockDiagonalSolver.h"
using namespace std;

enum class ChebyshevKind : unsigned
{
	// Chebyshev polynomial of the first kind, minimal on [lambdaMax / Ratio, lambdaMax]: the classical Chebyshev smoother
	// (Adams, Brezina, Hu, Tuminaro, J. Comput. Phys. 188, 2003). The low end of the spectrum is left to the coarse levels.
	First,
	// Chebyshev polynomial of the fourth kind, which minimizes max sqrt(lambda / lambdaMax) |p(lambda)| on [0, lambdaMax],
	// the smoothing factor of the convergence bound of the V-cycle: no lower bound to choose (Lottes, Numer. Linear
	// Algebra Appl. 30, 2023; Phillips, Fischer, arXiv:2210.03179). Degree 1: Jacobi with omega = 4 / (3 lambdaMax).
	Fourth
};

// Chebyshev smoother preconditioned by the block Jacobi iteration: after m iterations (MaxIterations: the degree of the
// polynomial), the error is p(D^-1 A) times the initial error, where D is the block diagonal of A and p the Chebyshev
// polynomial of the kind, with p(0) = 1. Each iteration costs a product by A and the solves of the diagonal blocks, as
// a block Jacobi iteration, but its block rows are independent: it runs in parallel, and its results don't depend on
// the number of threads. Each pass over the matrix fuses the product with the vector updates of the iteration.
// lambdaMax: the largest eigenvalue of D^-1 A, estimated in the setup by Lanczos, times EigenvalueSafetyFactor (the
// estimate is from below; above lambdaMax, |p| grows quickly).
// Requires a symmetric positive (semi-)definite matrix (real spectrum of D^-1 A); the smoother is then symmetric in the
// A-inner product: the same in pre- and post-smoothing.
class BlockChebyshev : public IterativeSolver
{
public:
	// Ratio lambdaMax / lambdaMin of the interval of the first kind (-cheb-ratio). 10: PETSc's interval, [0.1, 1.1] x the
	// estimate; better than Adams et al.'s 30 with U-AMG (PERFORMANCE.md)
	inline static double DefaultRatio = 10;
	// Lanczos iterations of the estimation of lambdaMax, and the factor applied to the estimate
	static constexpr int EigenvalueEstimateIterations = 10;
	static constexpr double EigenvalueSafetyFactor = 1.1;

private:
	int _blockSize;
	ChebyshevKind _kind;
	double _ratio;
	double _lambdaMaxEstimate = 0; // Lanczos estimate of the largest eigenvalue of D^-1 A (0 before the setup)

	BlockDiagonalSolver _diagBlockSolver; // blocks of size > 1
	Vector _invDiag;                      // blocks of size 1: the inverses of the diagonal coefficients
	Vector _d, _dNext, _r;                // work vectors of the iterations, kept from one call to the next
public:

	BlockChebyshev(int blockSize, ChebyshevKind kind, double ratio = DefaultRatio)
	{
		_blockSize = blockSize;
		_kind = kind;
		_ratio = ratio;
		if (_kind == ChebyshevKind::First && !(_ratio > 1))
			Utils::FatalError("Chebyshev smoother: the ratio lambdaMax / lambdaMin of its interval must be > 1, got " + to_string(_ratio) + ".");
	}

	virtual void Serialize(ostream& os) const override
	{
		os << (_blockSize == 1 ? "Chebyshev" : "Block Chebyshev") << " (";
		if (_blockSize != 1)
			os << "blockSize=" << _blockSize << ", ";
		if (_kind == ChebyshevKind::First)
			os << "1st kind on [lambdaMax/" << _ratio << ", lambdaMax]";
		else
			os << "4th kind on [0, lambdaMax]";
		if (_lambdaMaxEstimate > 0)
		{
			auto precision = os.precision(3);
			os << ", lambdaMax = " << EigenvalueSafetyFactor << " x " << _lambdaMaxEstimate;
			os.precision(precision);
		}
		os << ")";
	}

	// Lanczos estimate of the largest eigenvalue of D^-1 A (without the safety factor)
	double LambdaMaxEstimate() const
	{
		return _lambdaMaxEstimate;
	}

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);
		if (!A.isCompressed())
			Utils::FatalError("The Chebyshev smoothers require a compressed matrix (see Eigen's makeCompressed()).");

		BigNumber nb = A.rows() / _blockSize;
		if (_blockSize == 1)
		{
			_invDiag = Vector(A.rows());
			bool positiveDiagonal = true;
			#pragma omp parallel for reduction(&&:positiveDiagonal)
			for (BigNumber i = 0; i < A.rows(); i++)
			{
				double a_ii = 0;
				for (RowMajorSparseMatrix::InnerIterator it(A, i); it; ++it)
				{
					if (it.col() == i)
						a_ii = it.value();
				}
				positiveDiagonal = positiveDiagonal && a_ii > 0;
				_invDiag[i] = 1 / a_ii;
			}
			if (!positiveDiagonal)
				Utils::FatalError("The Chebyshev smoothers require a symmetric positive definite matrix: non-positive diagonal coefficient.");
		}
		else
			_diagBlockSolver.Setup(A, _blockSize);

		_lambdaMaxEstimate = EstimateLambdaMax();

		Flops lanczosIteration = Cost::MatVec(A) + 2 * _blockSize * A.rows() + 10 * A.rows();
		this->SetupComputationalWork = (nb * 2.0/3.0*pow(_blockSize, 3) + EigenvalueEstimateIterations * lanczosIteration) * 1e-6;
	}

	// The diagonal blocks and the estimate of lambdaMax of another Chebyshev smoother set up for A with the same block
	// size, whatever its kind
	bool SetupFrom(const SparseMatrix& A, const IterativeSolver& other) override
	{
		const BlockChebyshev* o = dynamic_cast<const BlockChebyshev*>(&other);
		if (!o || o->Matrix != &A || o->_blockSize != _blockSize || this->ComputeExactSolution)
			return false;
		IterativeSolver::Setup(A);
		_diagBlockSolver = o->_diagBlockSolver;
		_invDiag = o->_invDiag;
		_lambdaMaxEstimate = o->_lambdaMaxEstimate;
		this->SetupComputationalWork = o->SetupComputationalWork;
		return true;
	}

	bool CanOptimizeResidualComputation() override
	{
		return true;
	}

	// m = MaxIterations iterations (the degree of the polynomial), with
	//   r_0 = b - A x_0,  d_0 = c_0 D^-1 r_0,
	//   x_k = x_{k-1} + d_{k-1},  r_k = r_{k-1} - A d_{k-1},  d_k = a_k d_{k-1} + c_k D^-1 r_k  (k = 1..m-1),
	//   x_m = x_{m-1} + d_{m-1}  (and r_m = r_{m-1} - A d_{m-1} for the residual or Ax).
	// Each line is one pass over the block rows. The residual is that of the recurrence: b - A x up to rounding.
	void Solve(const Vector& b, Vector& x, bool xEquals0, bool computeResidual, bool computeAx) override
	{
		if (this->MaxIterations == 0)
		{
			IterativeSolver::Solve(b, x, xEquals0, computeResidual, computeAx);
			return;
		}
		assert(this->StoppingCrit == StoppingCriteria::MaxIterations); // a smoother only: a fixed polynomial

		const SparseMatrix& A = *this->Matrix;
		BigNumber n = A.rows();
		int m = this->MaxIterations;
		if (xEquals0 && x.rows() != n)
			x.resize(n);

		vector<double> a, c;
		Coefficients(m, a, c);

		_d.resize(n);
		_dNext.resize(n);
		_r.resize(n);
		this->Residual = Vector();
		Flops matVec = Cost::MatVec(A);
		Flops update = (2 * _blockSize + 3) * n; // D^-1 r and the linear combination of d
		Flops flops = 0;

		const double* r = b.data(); // r_k (b while x = 0)
		if (xEquals0)
		{
			Pass(nullptr, b.data(), nullptr, nullptr, _d.data(), 0, c[0], nullptr, false);                         flops += update;
		}
		else
		{
			Pass(x.data(), b.data(), _r.data(), nullptr, _d.data(), 0, c[0], nullptr, false);                      flops += matVec + update;
			r = _r.data();
		}
		for (int k = 1; k < m; k++)
		{
			Pass(_d.data(), r, _r.data(), _d.data(), _dNext.data(), a[k], c[k], x.data(), xEquals0);              flops += matVec + update + n;
			r = _r.data();
			xEquals0 = false;
			_d.swap(_dNext);
		}
		bool residual = computeResidual || computeAx;
		Pass(residual ? _d.data() : nullptr, r, _r.data(), _d.data(), nullptr, 0, 0, x.data(), xEquals0);         flops += (residual ? matVec : 0) + n;

		if (computeAx)
		{
			this->Ax = b - _r;                                                                                    flops += n;
		}
		if (computeResidual)
			this->Residual = std::move(_r);

		this->IterationCount = m;
		this->SolvingComputationalWork = flops * 1e-6;
	}

private:
	// The coefficients a_k (a_0 = 0, unused) and c_k of the m iterations, from lambdaMax = EigenvalueSafetyFactor x the
	// estimate
	void Coefficients(int m, vector<double>& a, vector<double>& c) const
	{
		double lambdaMax = EigenvalueSafetyFactor * _lambdaMaxEstimate;
		a.assign(m, 0);
		c.assign(m, 0);
		if (_kind == ChebyshevKind::First)
		{
			// Saad, Iterative Methods for Sparse Linear Systems, 2nd ed., 2003, Algorithm 12.1, on [lambdaMin, lambdaMax]:
			// error polynomial T_m((theta - lambda) / delta) / T_m(theta / delta)
			double lambdaMin = lambdaMax / _ratio;
			double theta = (lambdaMax + lambdaMin) / 2;
			double delta = (lambdaMax - lambdaMin) / 2;
			double sigma = theta / delta;
			double rho = 1 / sigma;
			c[0] = 1 / theta;
			for (int k = 1; k < m; k++)
			{
				double rhoNext = 1 / (2 * sigma - rho);
				a[k] = rhoNext * rho;
				c[k] = 2 * rhoNext / delta;
				rho = rhoNext;
			}
		}
		else
		{
			// Lottes 2023: error polynomial W_m(1 - 2 lambda / lambdaMax) / (2m + 1), W_m the Chebyshev polynomial of the
			// fourth kind
			c[0] = 4 / (3 * lambdaMax);
			for (int k = 1; k < m; k++)
			{
				a[k] = (2.0 * k - 1) / (2.0 * k + 3);
				c[k] = (8.0 * k + 4) / ((2.0 * k + 3) * lambdaMax);
			}
		}
	}

	// One pass over the block rows I, independent of each other:
	//   r_I = rIn_I - (A v)_I               if v (otherwise r = rIn, and r is not written),
	//   dOut_I = a dIn_I + c D_I^-1 r_I     if dOut (a dIn_I omitted if dIn is null),
	//   x_I = dIn_I (xEquals0) or x_I += dIn_I   if x.
	// rIn may be r (in place). v must be neither r nor dOut nor x.
	void Pass(const double* v, const double* rIn, double* r, const double* dIn, double* dOut, double a, double c, double* x, bool xEquals0)
	{
		const SparseMatrix& A = *this->Matrix;
		const SparseMatrixIndex* outer = A.outerIndexPtr();
		const SparseMatrixIndex* col = A.innerIndexPtr();
		const double* val = A.valuePtr();
		BigNumber nb = A.rows() / _blockSize;

		if (_blockSize == 1)
		{
			const double* invDiag = _invDiag.data();
			#pragma omp parallel for if (A.nonZeros() > 20000)
			for (BigNumber i = 0; i < nb; i++)
			{
				if (v || dOut)
				{
					double r_i = rIn[i];
					if (v)
					{
						for (SparseMatrixIndex p = outer[i]; p < outer[i + 1]; p++)
							r_i -= val[p] * v[col[p]];
						r[i] = r_i;
					}
					if (dOut)
						dOut[i] = (dIn ? a * dIn[i] : 0) + c * (invDiag[i] * r_i);
				}
				if (x)
					x[i] = xEquals0 ? dIn[i] : x[i] + dIn[i];
			}
		}
		else
		{
			#pragma omp parallel if (A.nonZeros() > 20000)
			{
				Vector r_I(_blockSize), z_I(_blockSize); // work vectors, allocated once per thread
				#pragma omp for
				for (BigNumber I = 0; I < nb; I++)
				{
					BigNumber i0 = I * _blockSize;
					for (int k = 0; k < _blockSize && (v || dOut); k++)
					{
						BigNumber i = i0 + k;
						double r_i = rIn[i];
						if (v)
						{
							for (SparseMatrixIndex p = outer[i]; p < outer[i + 1]; p++)
								r_i -= val[p] * v[col[p]];
							r[i] = r_i;
						}
						r_I[k] = r_i;
					}
					if (dOut)
					{
						_diagBlockSolver.Solve(I, r_I, z_I); // r_I destroyed
						for (int k = 0; k < _blockSize; k++)
							dOut[i0 + k] = (dIn ? a * dIn[i0 + k] : 0) + c * z_I[k];
					}
					if (x)
					{
						for (int k = 0; k < _blockSize; k++)
							x[i0 + k] = xEquals0 ? dIn[i0 + k] : x[i0 + k] + dIn[i0 + k];
					}
				}
			}
		}
	}

	// z = D^-1 r
	void ApplyInverseDiagonal(const Vector& r, Vector& z) const
	{
		const SparseMatrix& A = *this->Matrix;
		if (_blockSize == 1)
		{
			z = _invDiag.cwiseProduct(r);
			return;
		}
		BigNumber nb = A.rows() / _blockSize;
		#pragma omp parallel if (A.nonZeros() > 20000)
		{
			Vector r_I(_blockSize), z_I(_blockSize);
			#pragma omp for
			for (BigNumber I = 0; I < nb; I++)
			{
				r_I = r.segment(I * _blockSize, _blockSize);
				_diagBlockSolver.Solve(I, r_I, z_I);
				z.segment(I * _blockSize, _blockSize) = z_I;
			}
		}
	}

	// Largest eigenvalue of D^-1 A, by Lanczos: D^-1 A is similar to the symmetric D^-1/2 A D^-1/2, whose Lanczos
	// tridiagonal matrix is given by the coefficients of the conjugate gradient preconditioned by D (Saad, Iterative
	// Methods for Sparse Linear Systems, 2nd ed., 2003, §6.7.3). Its largest eigenvalue approaches lambdaMax from below.
	// Pseudo-random start vector, a function of the index only: the estimate is the same in every run, whatever the
	// number of threads.
	double EstimateLambdaMax() const
	{
		const SparseMatrix& A = *this->Matrix;
		BigNumber n = A.rows();
		Vector r(n), z(n), p, Ap(n);
		for (BigNumber i = 0; i < n; i++)
			r[i] = PseudoRandom(i);
		ApplyInverseDiagonal(r, z);
		double rz = r.dot(z);
		if (!(rz > 0))
			Utils::FatalError("The Chebyshev smoothers require a symmetric positive definite matrix: its diagonal blocks are not positive definite.");
		double rz0 = rz;
		p = z;
		vector<double> alpha, beta;
		for (int j = 0; j < EigenvalueEstimateIterations; j++)
		{
			Ap.noalias() = A * p;
			double pAp = p.dot(Ap);
			if (!(pAp > 0))
				Utils::FatalError("The Chebyshev smoothers require a symmetric positive definite matrix: p^T A p = " + to_string(pAp) + " in the estimation of its largest eigenvalue.");
			alpha.push_back(rz / pAp);
			if (j == EigenvalueEstimateIterations - 1)
				break;
			r -= alpha.back() * Ap;
			ApplyInverseDiagonal(r, z);
			double rzNext = r.dot(z);
			if (!(rzNext > 1e-28 * rz0)) // invariant Krylov subspace: the eigenvalues of the tridiagonal matrix are exact
				break;
			beta.push_back(rzNext / rz);
			p = z + beta.back() * p;
			rz = rzNext;
		}

		int m = alpha.size();
		DenseMatrix T = DenseMatrix::Zero(m, m);
		for (int j = 0; j < m; j++)
		{
			T(j, j) = 1 / alpha[j] + (j > 0 ? beta[j - 1] / alpha[j - 1] : 0);
			if (j > 0)
				T(j, j - 1) = T(j - 1, j) = sqrt(beta[j - 1]) / alpha[j - 1];
		}
		Eigen::SelfAdjointEigenSolver<DenseMatrix> eigenSolver(T, Eigen::EigenvaluesOnly);
		return eigenSolver.eigenvalues().maxCoeff();
	}

	// In [-1, 1), from the index (SplitMix64)
	static double PseudoRandom(uint64_t i)
	{
		uint64_t z = i + 0x9E3779B97F4A7C15ull;
		z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
		z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
		z ^= z >> 31;
		return (z >> 11) * 0x1.0p-52 - 1;
	}
};
