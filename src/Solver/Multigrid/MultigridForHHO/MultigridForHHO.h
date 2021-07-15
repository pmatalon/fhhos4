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
	GMG_H_Prolongation H_Prolongation = GMG_H_Prolongation::CellInterp_Trace;
	GMG_P_Prolongation P_Prolongation = GMG_P_Prolongation::Injection;
	GMG_P_Restriction P_Restriction = GMG_P_Restriction::RemoveHigherOrders;
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
		if (this->HP_CS == HP_CoarsStgy::HP_then_H)
			os << "hp then h [-hp-cs hp_h]";
		else if (this->HP_CS == HP_CoarsStgy::H_only)
			os << "h only [-hp-cs h]";
		else if (this->HP_CS == HP_CoarsStgy::P_only)
			os << "p only [-hp-cs p]";
		else if (this->HP_CS == HP_CoarsStgy::P_then_H)
			os << "p then h [-hp-cs p_h]";
		else if (this->HP_CS == HP_CoarsStgy::P_then_HP)
			os << "p then hp [-hp-cs p_hp]";
		os << endl;

		if (this->HP_CS == HP_CoarsStgy::P_only || this->HP_CS == HP_CoarsStgy::P_then_H || this->HP_CS == HP_CoarsStgy::P_then_HP)
		{
			os << "\t" << "p-prolongation          : ";
			if (P_Prolongation == GMG_P_Prolongation::Injection)
				os << "natural injection ";
			else if (P_Prolongation == GMG_P_Prolongation::H_Prolongation)
				os << "h-prolongation ";
			os << "[-p-prolong " << (unsigned)P_Prolongation << "]" << endl;

			os << "\t" << "p-restriction           : ";
			if (P_Restriction == GMG_P_Restriction::RemoveHigherOrders)
				os << "remove higher-order components ";
			else if (P_Restriction == GMG_P_Restriction::P_Transpose)
				os << "transpose of p-prolongation ";
			os << "[-p-restrict " << (unsigned)P_Restriction << "]" << endl;
		}
		if (this->HP_CS == HP_CoarsStgy::H_only || this->HP_CS == HP_CoarsStgy::HP_then_H || this->HP_CS == HP_CoarsStgy::P_then_H || this->HP_CS == HP_CoarsStgy::P_then_HP)
		{
			if (this->HP_CS == HP_CoarsStgy::H_only || this->HP_CS == HP_CoarsStgy::P_then_H)
				os << "\t" << "h-prolongation          : ";
			else
				os << "\t" << "hp-prolongation         : ";

			if (H_Prolongation == GMG_H_Prolongation::CellInterp_Trace)
				os << "coarse cell interpolation + trace/proj. on fine faces ";
			else if (H_Prolongation == GMG_H_Prolongation::CellInterp_Inject_Trace)
				os << "coarse cell interpolation + injection coarse to fine cells + trace on fine faces ";
			else if (H_Prolongation == GMG_H_Prolongation::CellInterp_ExactL2proj_Trace)
				os << "coarse cell interpolation + exact L2-proj. to fine cells + trace on fine faces ";
			else if (H_Prolongation == GMG_H_Prolongation::CellInterp_ApproxL2proj_Trace)
				os << "coarse cell interpolation + approx. L2-proj. to fine cells + trace on fine faces ";
			else if (H_Prolongation == GMG_H_Prolongation::CellInterp_FinerApproxL2proj_Trace)
				os << "coarse cell interpolation + approx. L2-proj. with subtriangulation + trace on fine faces ";
			else if (H_Prolongation == GMG_H_Prolongation::CellInterp_InjectAndTrace)
				os << "injection for common faces, and coarse cell interpolation + trace for the other ";
			else if (H_Prolongation == GMG_H_Prolongation::CellInterp_Inject_Adjoint)
				os << "coarse cell interpolation + injection coarse to fine cells + adjoint of cell interpolation ";
			else if (H_Prolongation == GMG_H_Prolongation::Wildey)
				os << "Wildey et al. ";
			else if (H_Prolongation == GMG_H_Prolongation::FaceInject)
				os << "injection coarse to fine faces ";
			os << "[-prolong " << (unsigned)H_Prolongation << "]" << endl;

			if (H_Prolongation != GMG_H_Prolongation::FaceInject && H_Prolongation != GMG_H_Prolongation::Wildey)
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
		if (Utils::RequiresNestedHierarchy(H_Prolongation) && !Utils::BuildsNestedMeshHierarchy(this->CoarseningStgy))
		{
			os << endl;
			Utils::Warning(os, "The selected coarsening strategy generates non-nested meshes, while the selected prolongation operator is made for nested meshes. Option -prolong " + to_string((unsigned)GMG_H_Prolongation::CellInterp_ExactL2proj_Trace) + ", " + to_string((unsigned)GMG_H_Prolongation::CellInterp_ApproxL2proj_Trace) + " or " + to_string((unsigned)GMG_H_Prolongation::CellInterp_FinerApproxL2proj_Trace) + " recommended.");
		}
		if (this->HP_CS == HP_CoarsStgy::P_only || this->HP_CS == HP_CoarsStgy::P_then_H || this->HP_CS == HP_CoarsStgy::P_then_HP)
		{
			if (P_Prolongation == GMG_P_Prolongation::Injection && (!this->_problem->HHO->FaceBasis->IsHierarchical || !this->_problem->HHO->OrthonormalizeBases))
				Utils::Warning("The natural injection for p-multigrid is implemented based on the assumption that the face bases are hierarchical and orthonormalized. Degraded convergence may be experienced.");
			if (P_Restriction == GMG_P_Restriction::RemoveHigherOrders && (!this->_problem->HHO->FaceBasis->IsHierarchical || !this->_problem->HHO->OrthonormalizeBases))
				Utils::Warning("The restriction for p-multigrid consisting in removing the higher-orders is implemented based on the assumption that the face bases are hierarchical and orthonormalized. Degraded convergence may be experienced.");
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
		return new LevelForHHO<Dim>(0, _problem, H_Prolongation, P_Prolongation, P_Restriction, UseHigherOrderReconstruction, UseHeterogeneousWeighting);
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
		
		LevelForHHO<Dim>* coarseLevel = new LevelForHHO<Dim>(fineLevel->Number + 1, coarseProblem, H_Prolongation, P_Prolongation, P_Restriction, UseHigherOrderReconstruction, UseHeterogeneousWeighting);
		return coarseLevel;
	}
};