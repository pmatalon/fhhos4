#pragma once
#include <mutex>
#include "UncondensedLevel.h"
#include "../Multigrid.h"
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
	Vector _cellConstants; // coordinate of the constant function 1 on the first basis function of each cell and face
	Vector _faceConstants; // (empty if they are all 1)
public:

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

	// cellConstants, faceConstants: coordinate of the constant function 1 on the first basis function of each cell and
	// face, the bases being hierarchical with a constant first function. Used to correct a restrictive assumption of
	// the paper, for which they are all 1 (see UncondensedLevel::CoarsenMesh()).
	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F, const Vector& cellConstants, const Vector& faceConstants) override
	{
		if (cellConstants.rows() * _cellBlockSize != A_T_F.rows() || faceConstants.rows() * _faceBlockSize != A_T_F.cols())
			Utils::FatalError("UncondensedAMG: one coordinate of the constant function is expected per cell (" + to_string(A_T_F.rows() / _cellBlockSize) + ") and per face (" + to_string(A_T_F.cols() / _faceBlockSize) + "), got " + to_string(cellConstants.rows()) + " and " + to_string(faceConstants.rows()) + ".");
		if (!cellConstants.allFinite() || !faceConstants.allFinite() || (cellConstants.array() == 0).any() || (faceConstants.array() == 0).any())
			Utils::FatalError("UncondensedAMG: the coordinates of the constant function must be finite and non-zero.");

		this->_fineLevel = this->CreateFineLevel();
		UncondensedLevel* fine = dynamic_cast<UncondensedLevel*>(this->_fineLevel);
		fine->A_T_T = &A_T_T;
		fine->A_T_F = &A_T_F;
		fine->A_F_F = &A_F_F;
		fine->LocalOp = &_localOperator;

		// With the coordinates of the paper (all 1), no rescaling: the setup is unchanged
		bool unitConstants = (cellConstants.array() == 1).all() && (faceConstants.array() == 1).all();
		_cellConstants = unitConstants ? Vector() : cellConstants;
		_faceConstants = unitConstants ? Vector() : faceConstants;
		fine->CellConstants = unitConstants ? nullptr : &_cellConstants;
		fine->FaceConstants = unitConstants ? nullptr : &_faceConstants;

		if (Utils::IsRefinementStrategy(this->H_CS))
			this->H_CS = H_CoarsStgy::MultiplePairwiseAggregation;
		if (this->H_CS == H_CoarsStgy::MultiplePairwiseAggregation && this->CoarseningFactor == 0)
			this->CoarseningFactor = 3.8;

		Multigrid::Setup(A);

		_localOperator.Current.Free();
		_localOperator.Next.Free();
		_localOperator.Assembled = nullptr;
	}

	Vector Solve(const Vector& b, string initialGuessCode) override
	{
		if (initialGuessCode.compare("smooth") == 0)
			Utils::Warning("Smooth initial guess unmanaged in AMG.");
		return Multigrid::Solve(b, initialGuessCode);
	}
	
private:
	// The coarse blocks A_F_F are not used by the algorithm of the paper (only A_T_T and A_T_F are), but by some options
	bool CoarseA_F_FNeeded() const
	{
		bool pLevelAfterHLevel = this->HP_CS == HP_CoarsStgy::H_then_P || this->HP_CS == HP_CoarsStgy::HP_then_P || this->HP_CS == HP_CoarsStgy::Alternate;
		return pLevelAfterHLevel // the p-coarsening extracts its blocks from those of the finer level
			|| this->FaceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsingAndTryAggregInteriorToInterfaces // face couplings
			|| _coarseningProlong == UAMGProlongation::HighOrder || _multigridProlong == UAMGProlongation::HighOrder; // trace on the removed faces
	}

	// Q_F chained over the coarsening steps is used by the multigrid prolongation if it is not the chained coarsening prolongations,
	// and by the coarse operator if it is not the Galerkin one
	bool ChainedQ_FNeeded() const
	{
		return _multigridProlong != UAMGProlongation::ChainedCoarseningProlongations || !this->UseGalerkinOperator;
	}

	UncondensedLevel* CreateLevel(int number, int degree, int cellBlockSize, int faceBlockSize) const
	{
		UncondensedLevel* level = new UncondensedLevel(number, degree, cellBlockSize, faceBlockSize, _strongCouplingThreshold, _faceProlong, _coarseningProlong, _multigridProlong);
		level->ComputeCoarseA_F_F = CoarseA_F_FNeeded();
		level->ComputeQ_F = ChainedQ_FNeeded();
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