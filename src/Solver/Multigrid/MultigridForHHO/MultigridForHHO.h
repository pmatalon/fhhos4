#pragma once
#include "../Multigrid.h"
#include "LevelForHHO.h"
using namespace std;

template <int Dim>
class MultigridForHHO : public Multigrid
{
private:
	Diffusion_HHO<Dim>* _problem;
public:
	GMGProlongation Prolongation = GMGProlongation::CellInterp_Trace;
	bool UseHigherOrderReconstruction = true;
	bool UseHeterogeneousWeighting = true;
	
	MultigridForHHO(Diffusion_HHO<Dim>* problem, int nLevels = 0)
		: Multigrid(MGType::h_Multigrid, nLevels)
	{
		this->InitializeWithProblem(problem);
	}

	// Use Initialize after these constructors
	MultigridForHHO() : MultigridForHHO(0)
	{}
	MultigridForHHO(int nLevels) : Multigrid(nLevels)
	{}

	void InitializeWithProblem(Diffusion_HHO<Dim>* problem)
	{
		this->_problem = problem;
		this->BlockSizeForBlockSmoothers = problem->HHO->nFaceUnknowns;

		if (problem->HHO->FaceBasis->GetDegree() == 0)
			Utils::Warning("k=0 is a special case for this multigrid method. Non-optimal convergence may be observed, especially in 3D.");
	}

	void BeginSerialize(ostream& os) const override
	{
		os << "MultigridForHHO" << endl;

		os << "\t" << "hp-strategy             : ";
		if (this->HP_Stgy == HP_Strategy::HP_then_H)
			os << "hp then h [-hp-stgy hp_h]";
		else if (this->HP_Stgy == HP_Strategy::H_only)
			os << "h only [-hp-stgy h]";
		else if (this->HP_Stgy == HP_Strategy::P_only)
			os << "p only [-hp-stgy p]";
		else if (this->HP_Stgy == HP_Strategy::P_then_H)
			os << "p then h [-hp-stgy p_h]";
		else if (this->HP_Stgy == HP_Strategy::P_then_HP)
			os << "p then hp [-hp-stgy p_hp]";
		os << endl;

		if (this->HP_Stgy == HP_Strategy::P_only || this->HP_Stgy == HP_Strategy::P_then_H || this->HP_Stgy == HP_Strategy::P_then_HP)
		{
			os << "\t" << "p-prolongation          : natural injection ";
			os << endl;
		}
		if (this->HP_Stgy == HP_Strategy::H_only || this->HP_Stgy == HP_Strategy::HP_then_H || this->HP_Stgy == HP_Strategy::P_then_H || this->HP_Stgy == HP_Strategy::P_then_HP)
		{
			if (this->HP_Stgy == HP_Strategy::H_only || this->HP_Stgy == HP_Strategy::P_then_H)
				os << "\t" << "h-prolongation          : ";
			else
				os << "\t" << "hp-prolongation         : ";

			if (Prolongation == GMGProlongation::CellInterp_Trace)
				os << "coarse cell interpolation + trace/proj. on fine faces ";
			else if (Prolongation == GMGProlongation::CellInterp_Inject_Trace)
				os << "coarse cell interpolation + injection coarse to fine cells + trace on fine faces ";
			else if (Prolongation == GMGProlongation::CellInterp_ExactL2proj_Trace)
				os << "coarse cell interpolation + exact L2-proj. to fine cells + trace on fine faces ";
			else if (Prolongation == GMGProlongation::CellInterp_ApproxL2proj_Trace)
				os << "coarse cell interpolation + approx. L2-proj. to fine cells + trace on fine faces ";
			else if (Prolongation == GMGProlongation::CellInterp_FinerApproxL2proj_Trace)
				os << "coarse cell interpolation + approx. L2-proj. with subtriangulation + trace on fine faces ";
			else if (Prolongation == GMGProlongation::CellInterp_InjectAndTrace)
				os << "injection for common faces, and coarse cell interpolation + trace for the other ";
			else if (Prolongation == GMGProlongation::CellInterp_Inject_Adjoint)
				os << "coarse cell interpolation + injection coarse to fine cells + adjoint of cell interpolation ";
			else if (Prolongation == GMGProlongation::Wildey)
				os << "Wildey et al. ";
			else if (Prolongation == GMGProlongation::FaceInject)
				os << "injection coarse to fine faces ";
			os << "[-prolong " << (unsigned)Prolongation << "]" << endl;

			if (Prolongation != GMGProlongation::FaceInject && Prolongation != GMGProlongation::Wildey)
			{
				if (!UseHigherOrderReconstruction)
					os << "\t" << "Higher-order reconstruc.: disabled" << endl;
				if (!UseHeterogeneousWeighting)
					os << "\t" << "Heterogeneous weighting : disabled" << endl;
			}
		}
	}

	void EndSerialize(ostream& os) const override
	{
		if (Utils::RequiresNestedHierarchy(Prolongation) && !Utils::BuildsNestedMeshHierarchy(this->CoarseningStgy))
		{
			os << endl;
			Utils::Warning(os, "The selected coarsening strategy generates non-nested meshes, while the selected prolongation operator is made for nested meshes. Option -prolong " + to_string((unsigned)GMGProlongation::CellInterp_ExactL2proj_Trace) + ", " + to_string((unsigned)GMGProlongation::CellInterp_ApproxL2proj_Trace) + " or " + to_string((unsigned)GMGProlongation::CellInterp_FinerApproxL2proj_Trace) + " recommended.");
		}
	}

	Vector Solve(const Vector& b, string initialGuessCode) override
	{
		if (initialGuessCode.compare("smooth") == 0)
		{
			Vector initialGuess;
			DomFunction coarseErrorFunction = [](const DomPoint& p)
			{
				double x = p.X;
				double y = p.Y;
				//return sin(4 * M_PI * x)*sin(4 * M_PI * y);
				//return x * (1 - x) * y*(1 - y);
				return sin(M_PI * x)*sin(M_PI * y);
			};

			Level* coarseLevel = this->_fineLevel;
			while (coarseLevel->CoarserLevel)
				coarseLevel = coarseLevel->CoarserLevel;

			Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(coarseLevel)->_problem;

			Vector coarseError = coarsePb->ProjectOnFaceDiscreteSpace(coarseErrorFunction);

			Level* fineLevel = coarseLevel;
			Vector prolongatedError = coarseError;
			while (fineLevel->FinerLevel)
			{
				fineLevel = fineLevel->FinerLevel;
				prolongatedError = fineLevel->Prolong(prolongatedError);
			}

			initialGuess = -prolongatedError;
			Multigrid::Solve(b, initialGuess, false);
			return initialGuess;
		}
		else
			return Multigrid::Solve(b, initialGuessCode);
	}
	
protected:
	Level* CreateFineLevel() const override
	{
		assert(_problem && "The multigrid has not been initialized with a problem");
		return new LevelForHHO<Dim>(0, _problem, Prolongation, UseHigherOrderReconstruction, UseHeterogeneousWeighting);
	}

	Level* CreateCoarseLevel(Level* fineLevel, CoarseningType coarseningType) override
	{
		LevelForHHO<Dim>* hhoFineLevel = dynamic_cast<LevelForHHO<Dim>*>(fineLevel);

		Diffusion_HHO<Dim>* coarseProblem;
		if (coarseningType == CoarseningType::H)
			coarseProblem = hhoFineLevel->_problem->GetProblemOnCoarserMesh();
		else if (coarseningType == CoarseningType::P)
			coarseProblem = hhoFineLevel->_problem->GetProblemForLowerDegree();
		else if (coarseningType == CoarseningType::HP)
			coarseProblem = hhoFineLevel->_problem->GetProblemOnCoarserMeshAndLowerDegree();
		else
			assert(false);
		
		LevelForHHO<Dim>* coarseLevel = new LevelForHHO<Dim>(fineLevel->Number + 1, coarseProblem, Prolongation, UseHigherOrderReconstruction, UseHeterogeneousWeighting);
		return coarseLevel;
	}
};