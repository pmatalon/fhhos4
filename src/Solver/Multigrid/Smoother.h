#pragma once
#include "../FixedPoint/BlockSOR.h"
#include "../FixedPoint/BlockGaussSeidel.h"
#include "../FixedPoint/BlockJacobi.h"
#include "../FixedPoint/GaussSeidel.h"
#include "../FixedPoint/BlockChebyshev.h"
using namespace std;

class Smoother
{
protected:
	IterativeSolver* _solver;
	int _nSmoothingIterations;
	FloatVector _singleResidual, _singleAx; // the residual and A x of the last smoothing in single precision

public:

	Smoother(IterativeSolver* solver, int nSmoothingIterations) :
		_solver(solver)
	{
		_solver->StoppingCrit = StoppingCriteria::MaxIterations;
		_solver->MaxIterations = nSmoothingIterations;
		_solver->PrintIterationResults = false;
		_solver->ComputeExactSolution = false;
	}

	int Iterations()
	{
		return _solver->MaxIterations;
	}

	void Setup(const SparseMatrix& A)
	{
		if (_solver->MaxIterations > 0)
			_solver->Setup(A);
	}

	// Same, reusing the setup of another smoother of A when it is of the same kind (the pre-smoother of the level)
	void Setup(const SparseMatrix& A, const Smoother& setUpForA)
	{
		if (_solver->MaxIterations > 0 && !(setUpForA._solver->MaxIterations > 0 && _solver->SetupFrom(A, *setUpForA._solver)))
			_solver->Setup(A);
	}

	// Mixed precision: with the matrix in single precision, for the smoothing of float vectors (see
	// SupportsSinglePrecision())
	void Setup(const FloatSparseMatrix& A)
	{
		if (_solver->MaxIterations > 0)
			_solver->SetupSinglePrecision(A);
	}

	void Setup(const FloatSparseMatrix& A, const Smoother& setUpForA)
	{
		if (_solver->MaxIterations > 0 && !(setUpForA._solver->MaxIterations > 0 && _solver->SetupSinglePrecisionFrom(A, *setUpForA._solver)))
			_solver->SetupSinglePrecision(A);
	}

	// On several processes
	void Setup(const DistributedMatrix& A)
	{
		if (_solver->MaxIterations > 0)
			_solver->Setup(A);
	}

	void Setup(const DistributedMatrix& A, const Smoother& setUpForA)
	{
		if (_solver->MaxIterations > 0 && !(setUpForA._solver->MaxIterations > 0 && _solver->SetupFrom(A, *setUpForA._solver)))
			_solver->Setup(A);
	}

	// The order of the (block) rows in all the sweeps of the Gauss-Seidel smoothers (see DownwindOrder()); the other
	// smoothers keep theirs
	void SetSweepOrder(const vector<BigNumber>* order)
	{
		if (GaussSeidel* gs = dynamic_cast<GaussSeidel*>(_solver))
			gs->SweepOrder = order;
		else if (BlockGaussSeidel* bgs = dynamic_cast<BlockGaussSeidel*>(_solver))
			bgs->SweepOrder = order;
	}

	// Mixed precision: the smoother in single precision (a smoother without iterations needs nothing)
	bool SupportsSinglePrecision() const
	{
		return _solver->MaxIterations == 0 || _solver->SupportsSinglePrecision();
	}

	bool CanOptimizeResidualComputation()
	{
		return _solver->MaxIterations > 0 && _solver->CanOptimizeResidualComputation(); // without iterations, no setup (the hybrid Gauss-Seidel's chunks need the matrix)
	}

	void Smooth(Vector& x, const Vector& b, bool& xEquals0)
	{
		_solver->Solve(b, x, xEquals0, false, false);
		xEquals0 = false;
	}

	void SmoothAndComputeResidualOrAx(Vector& x, const Vector& b, bool& xEquals0, bool computeResidual, bool computeAx)
	{
		_solver->Solve(b, x, xEquals0, computeResidual, computeAx);
		xEquals0 = false;
	}

	// The same in single precision
	void Smooth(FloatVector& x, const FloatVector& b, bool& xEquals0)
	{
		if (_solver->MaxIterations > 0)
			_solver->SolveInSinglePrecision(b, x, xEquals0, nullptr, nullptr);
		else if (xEquals0)
			x = FloatVector::Zero(b.rows());
		xEquals0 = false;
	}

	void SmoothAndComputeResidualOrAx(FloatVector& x, const FloatVector& b, bool& xEquals0, bool computeResidual, bool computeAx)
	{
		_solver->SolveInSinglePrecision(b, x, xEquals0, computeResidual ? &_singleResidual : nullptr, computeAx ? &_singleAx : nullptr);
		xEquals0 = false;
	}

	// The residual and A x of the last smoothing, in precision T
	template <typename T>
	VectorOf<T>&& TakeResidual()
	{
		if constexpr (std::is_same_v<T, float>)
			return std::move(_singleResidual);
		else
			return Residual();
	}

	template <typename T>
	VectorOf<T>&& TakeAx()
	{
		if constexpr (std::is_same_v<T, float>)
			return std::move(_singleAx);
		else
			return Ax();
	}

	friend ostream& operator<<(ostream& os, const Smoother& s)
	{
		os << *(s._solver);
		return os;
	}

	MFlops SolvingComputationalWork()
	{
		return _solver->SolvingComputationalWork;
	}

	Vector&& Residual()
	{
		return std::move(_solver->Residual);
	}

	Vector&& Ax()
	{
		return std::move(_solver->Ax);
	}

	virtual ~Smoother()
	{
		delete _solver;
	}

};

class BlockJacobiSmoother : public Smoother
{
public:
	BlockJacobiSmoother(int blockSize, double omega, int nSmoothingIterations) : Smoother(new BlockJacobi(blockSize, omega), nSmoothingIterations) {}
};

class BlockSORSmoother : public Smoother
{
public:
	BlockSORSmoother(int blockSize, double omega, Direction direction, int nSmoothingIterations) : Smoother(new BlockSOR(blockSize, omega, direction), nSmoothingIterations) {}
};

class BlockGaussSeidelSmoother : public Smoother
{
public:
	BlockGaussSeidelSmoother(int blockSize, Direction direction, bool hybrid, int nSmoothingIterations) : Smoother(new BlockGaussSeidel(blockSize, direction, hybrid), nSmoothingIterations) {}
};

// Number of smoothing iterations: the degree of the polynomial
class ChebyshevSmoother : public Smoother
{
public:
	ChebyshevSmoother(int blockSize, ChebyshevKind kind, int nSmoothingIterations) : Smoother(new BlockChebyshev(blockSize, kind), nSmoothingIterations) {}
};

class GaussSeidelSmoother : public Smoother
{
public:
	GaussSeidelSmoother(Direction direction, int nSmoothingIterations) : Smoother(new GaussSeidel(direction), nSmoothingIterations) {}
	GaussSeidelSmoother(Direction direction, bool hybrid, int nSmoothingIterations) : Smoother(new GaussSeidel(direction, hybrid), nSmoothingIterations) {}
};