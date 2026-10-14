#pragma once
#include <mutex>
#include "UncondensedLevel.h"
#include "../Multigrid.h"
#include "../../Direct/EigenSparseCholesky.h"
using namespace std;

class UncondensedAMG : public Multigrid
{
private:
	UAMGFaceProlongation _faceProlong = UAMGFaceProlongation::FaceAggregates;
	UAMGProlongation _coarseningProlong = UAMGProlongation::FaceProlongation;
	UAMGProlongation _multigridProlong = UAMGProlongation::ReconstructSmoothedTraceOrInject;
	int _dim;
	int _degree;
	int _cellBlockSize;
	int _faceBlockSize;
	double _strongCouplingThreshold;
	LocalOperator _localOperator; // local matrices of the operators of the coarsening passes, during the setup
	bool _symmetric = true;       // set by Setup(): see UncondensedLevel::Symmetric
public:
	bool ManageAnisotropy = true; // see UncondensedLevel::ManageAnisotropy
	bool DownwindOrdering = true; // see UncondensedLevel::DownwindOrdering


	UncondensedAMG(int dim, int degree, int cellBlockSize, int faceBlockSize, double strongCouplingThreshold, UAMGFaceProlongation faceProlong, UAMGProlongation coarseningProlong, UAMGProlongation mgProlong, int nLevels = 0)
		: Multigrid(nLevels)
	{
		this->_dim = dim;
		this->_degree = degree;
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
		this->_strongCouplingThreshold = strongCouplingThreshold;
		this->_faceProlong = faceProlong;
		this->_coarseningProlong = coarseningProlong;
		this->_multigridProlong = mgProlong;
		this->BlockSizeForBlockSmoothers = faceBlockSize;
		this->UseGalerkinOperator = true;
		this->H_CS = H_CoarsStgy::MultiplePairwiseAggregation;
		this->FaceCoarseningStgy = FaceCoarseningStrategy::InterfaceCollapsing;
	}

	void BeginSerialize(ostream& os) const override
	{
		os << "UncondensedAMG" << endl;

		os << "\t" << "Face prolongation       : ";
		if (_faceProlong == UAMGFaceProlongation::BoundaryAggregatesInteriorAverage)
			os << "aggregate interface faces, avg on interior faces ";
		else if (_faceProlong == UAMGFaceProlongation::BoundaryAggregatesInteriorZero)
			os << "aggregate interface faces, 0 on interior faces ";
		else if (_faceProlong == UAMGFaceProlongation::FaceAggregates)
			os << "aggregate all faces ";
		os << "[-face-prolong " << (unsigned)_faceProlong << "]" << endl;

		os << "\t" << "Coarsening Prolongation : ";
		if (_coarseningProlong == UAMGProlongation::ReconstructionTrace)
			os << "ReconstructionTrace ";
		else if (_coarseningProlong == UAMGProlongation::FaceProlongation)
			os << "FaceProlongation ";
		else if (_coarseningProlong == UAMGProlongation::ReconstructTraceOrInject)
			os << "ReconstructTraceOrInject ";
		else if (_coarseningProlong == UAMGProlongation::ReconstructSmoothedTraceOrInject)
			os << "ReconstructSmoothedTraceOrInject ";
		os << "[-coarsening-prolong " << (unsigned)_coarseningProlong << "]" << endl;

		os << "\t" << "Multigrid prolongation  : ";
		if (_multigridProlong == UAMGProlongation::ReconstructionTrace)
			os << "ReconstructionTrace ";
		else if (_multigridProlong == UAMGProlongation::ChainedCoarseningProlongations)
			os << "Chained coarsening prolongations ";
		else if (_multigridProlong == UAMGProlongation::FaceProlongation)
			os << "FaceProlongation ";
		else if (_multigridProlong == UAMGProlongation::ReconstructTraceOrInject)
			os << "ReconstructTraceOrInject ";
		else if (_multigridProlong == UAMGProlongation::ReconstructSmoothedTraceOrInject)
			os << "ReconstructSmoothedTraceOrInject ";
		os << "[-prolong " << (unsigned)_multigridProlong << "]" << endl;
	}

	void EndSerialize(ostream& os) const override
	{
	}

	void Setup(const SparseMatrix& A) override
	{
		Utils::FatalError("The method Setup(const SparseMatrix& A) cannot be used for this solver.");
	}

	// A_F_F: may be empty (0 x 0) if A_F_FNeeded() is false: the algorithm of the paper only uses A_T_T and A_T_F.
	// A_F_T: empty (0 x 0) if the matrix is symmetric (A_F_T = A_T_F^T). Otherwise, the same algorithm runs on the
	// matrix itself, all the blocks read in full (see UncondensedLevel::Symmetric): the couplings of the aggregation
	// from the averaged traces of A_T_F and A_F_T, the reconstruction -A_T_T^{-1} A_T_F (inverse of the static
	// condensation), the smoothing of the prolongation and the Galerkin coarse operators P^T A P with the matrix; then
	// the Gauss-Seidel smoothers follow the downwind order of each level (see DownwindOrder()), the K-cycle must iterate
	// with GCR, and the coarse solver must not be Cholesky.
	// cellInterpOfOne, faceInterpOfOne: interpolation of the function 1 on the polynomial bases of the cells and faces,
	// one coefficient c per cell and per face (numbered as the blocks of rows and columns of A_T_F). The bases must be
	// hierarchical with a constant first function phi_0: 1 = c phi_0. Used to correct a restrictive assumption of the
	// paper, where c = 1 everywhere (see UncondensedLevel::CoarsenMesh()).
	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_T, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
	{
		if (cellInterpOfOne.rows() != A_T_F.rows() / _cellBlockSize || faceInterpOfOne.rows() != A_T_F.cols() / _faceBlockSize)
			Utils::FatalError("UncondensedAMG: the interpolations of 1 on the cells and faces must have one coefficient per cell (" + to_string(A_T_F.rows() / _cellBlockSize) + ") and per face (" + to_string(A_T_F.cols() / _faceBlockSize) + "), got " + to_string(cellInterpOfOne.rows()) + " and " + to_string(faceInterpOfOne.rows()) + ".");
		if (!cellInterpOfOne.allFinite() || (cellInterpOfOne.array() == 0).any() || !faceInterpOfOne.allFinite() || (faceInterpOfOne.array() == 0).any())
			Utils::FatalError("UncondensedAMG: the interpolations of 1 on the cells and faces must be finite and non-zero.");
		_symmetric = A_F_T.rows() == 0;
		if (!_symmetric)
		{
			if (A_F_T.rows() != A_T_F.cols() || A_F_T.cols() != A_T_F.rows())
				Utils::FatalError("UncondensedAMG: A_F_T must have the sizes of the transpose of A_T_F.");
			if (this->Cycle == 'K' && this->KCycleKrylov != KrylovMethod::GCR)
				Utils::FatalError("UncondensedAMG: the K-cycle of a non-symmetric matrix must iterate with GCR, not FCG.");
			if (dynamic_cast<EigenSparseCholesky*>(this->CoarseSolver))
				Utils::FatalError("UncondensedAMG: the Cholesky factorization cannot be the coarse solver of a non-symmetric matrix (use LU).");
		}

		this->_fineLevel = this->CreateFineLevel();
		UncondensedLevel* fine = dynamic_cast<UncondensedLevel*>(this->_fineLevel);
		fine->A_T_T = &A_T_T;
		fine->A_T_F = &A_T_F;
		fine->A_F_T = _symmetric ? nullptr : &A_F_T;
		if (A_F_F.rows() == 0 && A_F_FNeeded())
			Utils::FatalError("UncondensedAMG: the chosen options need the block A_F_F.");
		fine->A_F_F = A_F_F.rows() > 0 ? &A_F_F : nullptr;
		fine->LocalOp = &_localOperator;

		// If the constant has coordinate 1 everywhere, as in the paper, no rescaling: the setup is unchanged
		if ((cellInterpOfOne.array() != 1).any() || (faceInterpOfOne.array() != 1).any())
		{
			fine->CellInterpOfOne = cellInterpOfOne;
			fine->FaceInterpOfOne = faceInterpOfOne;
		}

		if (Utils::IsRefinementStrategy(this->H_CS))
			this->H_CS = H_CoarsStgy::MultiplePairwiseAggregation;
		if (this->H_CS == H_CoarsStgy::MultiplePairwiseAggregation && this->CoarseningFactor == 0)
			this->CoarseningFactor = 3.8;

		Multigrid::Setup(A);

		_localOperator.Current.Free();
		_localOperator.Next.Free();
		_localOperator.Assembled = nullptr;

		// The blocks are only read during the setup: the caller may free them (the library fhhos4 does)
		fine->A_T_T = nullptr;
		fine->A_T_F = nullptr;
		fine->A_F_F = nullptr;
		fine->A_F_T = nullptr;
	}

	Vector Solve(const Vector& b, string initialGuessCode) override
	{
		if (initialGuessCode.compare("smooth") == 0)
			Utils::Warning("Smooth initial guess unmanaged in AMG.");
		return Multigrid::Solve(b, initialGuessCode);
	}
	
	// The blocks A_F_F are not used by the algorithm of the paper (only A_T_T and A_T_F are), but by some options
	bool A_F_FNeeded() const
	{
		bool pLevelAfterHLevel = this->HP_CS == HP_CoarsStgy::H_then_P || this->HP_CS == HP_CoarsStgy::HP_then_P || this->HP_CS == HP_CoarsStgy::Alternate;
		return pLevelAfterHLevel // the p-coarsening extracts its blocks from those of the finer level
			|| this->FaceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsingAndTryAggregInteriorToInterfaces // face couplings
			|| _coarseningProlong == UAMGProlongation::HighOrder || _multigridProlong == UAMGProlongation::HighOrder; // trace on the removed faces
	}

private:
	// Q_F chained over the coarsening steps is used by the multigrid prolongation if it is not the chained coarsening prolongations,
	// and by the coarse operator if it is not the Galerkin one
	bool ChainedQ_FNeeded() const
	{
		return _multigridProlong != UAMGProlongation::ChainedCoarseningProlongations || !this->UseGalerkinOperator;
	}

	UncondensedLevel* CreateLevel(int number, int degree, int cellBlockSize, int faceBlockSize) const
	{
		UncondensedLevel* level = new UncondensedLevel(number, degree, cellBlockSize, faceBlockSize, _strongCouplingThreshold, _faceProlong, _coarseningProlong, _multigridProlong);
		level->ComputeCoarseA_F_F = A_F_FNeeded();
		level->ComputeQ_F = ChainedQ_FNeeded();
		level->ManageAnisotropy = ManageAnisotropy;
		level->Symmetric = _symmetric;
		level->DownwindOrdering = DownwindOrdering;
		return level;
	}

protected:
	Level* CreateFineLevel() const override
	{
		return CreateLevel(0, _degree, _cellBlockSize, _faceBlockSize);
	}

	Level* CreateCoarseLevel(Level* fineLevel, CoarseningType coarseningType, int coarseDegree) override
	{
		if (coarseningType == CoarseningType::HP)
			Utils::FatalError("hp-coarsening not allowed for this multigrid.");

		UncondensedLevel* fine = dynamic_cast<UncondensedLevel*>(fineLevel);

		UncondensedLevel* coarse;
		if (coarseningType == CoarseningType::P)
		{
			int coarseCellBlockSize = Utils::Binomial(coarseDegree + _dim    , coarseDegree);
			int coarseFaceBlockSize = Utils::Binomial(coarseDegree + _dim - 1, coarseDegree);
			coarse = CreateLevel(fine->Number + 1, coarseDegree, coarseCellBlockSize, coarseFaceBlockSize);
		}
		else
			coarse = CreateLevel(fine->Number + 1, fine->PolynomialDegree(), fine->CellBlockSize(), fine->FaceBlockSize());
		coarse->LocalOp = fine->LocalOp;
		return coarse;
	}
};