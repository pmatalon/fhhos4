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
	int _cellBlockSize;
	int _faceBlockSize;
	double _strongCouplingThreshold;
public:

	UncondensedAMG(int cellBlockSize, int faceBlockSize, double strongCouplingThreshold, UAMGFaceProlongation faceProlong, UAMGProlongation coarseningProlong, UAMGProlongation mgProlong, int nLevels = 0)
		: Multigrid(nLevels)
	{
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
		assert(false && "This Setup method cannot be used in this solver.");
	}

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F) override
	{
		this->_fineLevel = this->CreateFineLevel();
		UncondensedLevel* fine = dynamic_cast<UncondensedLevel*>(this->_fineLevel);
		fine->A_T_T = &A_T_T;
		fine->A_T_F = &A_T_F;
		fine->A_F_F = &A_F_F;
		SparseMatrix* inv_A_T_T = new SparseMatrix(Utils::InvertBlockDiagMatrix(A_T_T, _cellBlockSize));
		fine->inv_A_T_T = inv_A_T_T;

		if (Utils::IsRefinementStrategy(this->H_CS))
			this->H_CS = H_CoarsStgy::MultiplePairwiseAggregation;

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
		return new UncondensedLevel(0, _cellBlockSize, _faceBlockSize, _strongCouplingThreshold, _faceProlong, _coarseningProlong, _multigridProlong);
	}

	Level* CreateCoarseLevel(Level* fineLevel, CoarseningType coarseningType) override
	{
		if (coarseningType != CoarseningType::H)
			Utils::FatalError("Only h-coarsening allowed for this multigrid.");

		UncondensedLevel* coarseLevel = new UncondensedLevel(fineLevel->Number + 1, _cellBlockSize, _faceBlockSize, _strongCouplingThreshold, _faceProlong, _coarseningProlong, _multigridProlong);
		coarseLevel->OperatorMatrix = &dynamic_cast<UncondensedLevel*>(fineLevel)->Ac;
		return coarseLevel;
	}
};