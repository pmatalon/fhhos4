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
public:
	bool ManageAnisotropy = true; // see UncondensedLevel::ManageAnisotropy


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
	// cellInterpOfOne, faceInterpOfOne: interpolation of the function 1 on the polynomial bases of the cells and faces,
	// numbered as the rows and columns of A_T_F. The bases must be hierarchical with a constant first function: the
	// interpolation of 1 is then c on the first DoF of each block, 0 on the others. Used to correct a restrictive
	// assumption of the paper, where c = 1 everywhere (see UncondensedLevel::CoarsenMesh()).
	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F, const Vector& cellInterpOfOne, const Vector& faceInterpOfOne) override
	{
		if (cellInterpOfOne.rows() != A_T_F.rows() || faceInterpOfOne.rows() != A_T_F.cols())
			Utils::FatalError("UncondensedAMG: the interpolations of 1 on the cells and faces must have the sizes of the rows (" + to_string(A_T_F.rows()) + ") and columns (" + to_string(A_T_F.cols()) + ") of A_T_F, got " + to_string(cellInterpOfOne.rows()) + " and " + to_string(faceInterpOfOne.rows()) + ".");
		CheckInterpOfOne(cellInterpOfOne, _cellBlockSize, "cell");
		CheckInterpOfOne(faceInterpOfOne, _faceBlockSize, "face");

		this->_fineLevel = this->CreateFineLevel();
		UncondensedLevel* fine = dynamic_cast<UncondensedLevel*>(this->_fineLevel);
		fine->A_T_T = &A_T_T;
		fine->A_T_F = &A_T_F;
		if (A_F_F.rows() == 0 && A_F_FNeeded())
			Utils::FatalError("UncondensedAMG: the chosen options need the block A_F_F.");
		fine->A_F_F = A_F_F.rows() > 0 ? &A_F_F : nullptr;
		fine->LocalOp = &_localOperator;

		// If the constant has coordinate 1 everywhere, as in the paper, no rescaling: the setup is unchanged
		if (!FirstCoeffsAreOne(cellInterpOfOne, _cellBlockSize) || !FirstCoeffsAreOne(faceInterpOfOne, _faceBlockSize))
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

		// The blocks are only read during the setup: the caller may free them (the library fhhos4_AMG does)
		fine->A_T_T = nullptr;
		fine->A_T_F = nullptr;
		fine->A_F_F = nullptr;
	}

	Vector Solve(const Vector& b, string initialGuessCode) override
	{
		if (initialGuessCode.compare("smooth") == 0)
			Utils::Warning("Smooth initial guess unmanaged in AMG.");
		return Multigrid::Solve(b, initialGuessCode);
	}
	
private:
	// The first coefficient of each block (the coordinate of the constant) must be finite and non-zero, the others
	// zero (hierarchical basis with the constant first)
	static void CheckInterpOfOne(const Vector& interpOfOne, int blockSize, string entity)
	{
		bool constantFirst = true;
		for (BigNumber i = 0; i < (BigNumber)interpOfOne.rows(); i += blockSize)
		{
			double c = interpOfOne[i];
			if (!std::isfinite(c) || c == 0)
				Utils::FatalError("UncondensedAMG: the interpolation of 1 on the " + entity + " bases must have a finite and non-zero first coefficient in each block.");
			for (int j = 1; j < blockSize; j++)
			{
				if (!(abs(interpOfOne[i + j]) <= 1e-6 * abs(c)))
					constantFirst = false;
			}
		}
		if (!constantFirst)
			Utils::Warning("UncondensedAMG: the interpolation of 1 on the " + entity + " bases is not zero beyond the first coefficient of each block: the bases are not hierarchical with the constant first, the coarse levels may not represent the constant functions.");
	}

	static bool FirstCoeffsAreOne(const Vector& interpOfOne, int blockSize)
	{
		for (BigNumber i = 0; i < (BigNumber)interpOfOne.rows(); i += blockSize)
		{
			if (interpOfOne[i] != 1)
				return false;
		}
		return true;
	}

public:
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