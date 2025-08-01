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

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F) override
	{
		this->_fineLevel = this->CreateFineLevel();
		UncondensedLevel* fine = dynamic_cast<UncondensedLevel*>(this->_fineLevel);
		fine->A_T_T = &A_T_T;
		fine->A_T_F = &A_T_F;
		fine->A_F_F = &A_F_F;

		if (Utils::IsRefinementStrategy(this->H_CS))
			this->H_CS = H_CoarsStgy::MultiplePairwiseAggregation;
		if (this->H_CS == H_CoarsStgy::MultiplePairwiseAggregation && this->CoarseningFactor == 0)
			this->CoarseningFactor = 3.8;

		Multigrid::Setup(A);
	}

	Vector Solve(const Vector& b, string initialGuessCode) override
	{
		if (initialGuessCode.compare("smooth") == 0)
			Utils::Warning("Smooth initial guess unmanaged in AMG.");
		return Multigrid::Solve(b, initialGuessCode);
	}
	
protected:
	Level* CreateFineLevel() const override
	{
		return new UncondensedLevel(0, _degree, _cellBlockSize, _faceBlockSize, _strongCouplingThreshold, _faceProlong, _coarseningProlong, _multigridProlong);
	}

	Level* CreateCoarseLevel(Level* fineLevel, CoarseningType coarseningType, int coarseDegree) override
	{
		if (coarseningType == CoarseningType::HP)
			Utils::FatalError("hp-coarsening not allowed for this multigrid.");

		UncondensedLevel* fine = dynamic_cast<UncondensedLevel*>(fineLevel);

		if (coarseningType == CoarseningType::P)
		{
			int coarseCellBlockSize = Utils::Binomial(coarseDegree + _dim    , coarseDegree);
			int coarseFaceBlockSize = Utils::Binomial(coarseDegree + _dim - 1, coarseDegree);
			UncondensedLevel* coarseLevel = new UncondensedLevel(fine->Number + 1, coarseDegree, coarseCellBlockSize, coarseFaceBlockSize, _strongCouplingThreshold, _faceProlong, _coarseningProlong, _multigridProlong);
			return coarseLevel;
		}
		else
		{
			UncondensedLevel* coarseLevel = new UncondensedLevel(fine->Number + 1, fine->PolynomialDegree(), fine->CellBlockSize(), fine->FaceBlockSize(), _strongCouplingThreshold, _faceProlong, _coarseningProlong, _multigridProlong);
			return coarseLevel;
		}
	}
};