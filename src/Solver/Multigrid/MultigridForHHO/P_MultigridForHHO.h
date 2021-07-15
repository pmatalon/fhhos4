#pragma once
#include "MultigridForHHO.h"
using namespace std;

template <int Dim>
class P_LevelForHHO : public Level
{
private:
public:
	Diffusion_HHO<Dim>* _problem;

	P_LevelForHHO(int number, Diffusion_HHO<Dim>* problem)
		: Level(number)
	{
		this->_problem = problem;
	}

	int PolynomialDegree() override
	{
		return _problem->HHO->FaceBasis->GetDegree();
	}

	int BlockSizeForBlockSmoothers() override
	{
		return _problem->HHO->nFaceUnknowns;
	}

	void ExportVector(const Vector& v, string suffix, int levelNumber) override
	{
		if (!this->IsFinestLevel())
			this->FinerLevel->ExportVector(v, suffix, levelNumber);
		else
			this->_problem->ExportVector(v, "level" + to_string(levelNumber) + "_" + suffix);
	}

	void ExportMatrix(const SparseMatrix& M, string suffix, int levelNumber) override
	{
		if (!this->IsFinestLevel())
			this->FinerLevel->ExportMatrix(M, suffix, levelNumber);
		else
			this->_problem->ExportMatrix(M, "level" + to_string(levelNumber) + "_" + suffix);
	}

	void SetupDiscretizedOperator() override
	{
		this->OperatorMatrix = &this->_problem->A;
	}

	void OnStartSetup() override
	{
		cout << "\t\tk = " << this->PolynomialDegree() << endl;

		if (!IsCoarsestLevel())
		{
			Diffusion_HHO<Dim>* coarsePb = dynamic_cast<P_LevelForHHO<Dim>*>(CoarserLevel)->_problem;

			if (this->UseGalerkinOperator)
				coarsePb->InitHHO();
			else
			{
				ActionsArguments actions;
				actions.LogAssembly = false;
				actions.AssembleRightHandSide = false;

				// If the bases are not hierarchical, then the bases for every degree have to be stored in ReferenceCartesianShape,
				// on the model of what is done in ReferenceShape with the use of maps.
				if (!_problem->HHO->FaceBasis->IsHierarchical)
					Utils::FatalError("This p-Multigrid is only implemented for the use of hierarchical bases.");
				//actions.InitReferenceShapes = this->MultigridType == MGType::p_Multigrid;

				coarsePb->Assemble(actions);
			}
		}
	}

	void SetupProlongation() override
	{
		Mesh<Dim>* mesh = this->_problem->_mesh;
		Diffusion_HHO<Dim>* higherPb = this->_problem;
		Diffusion_HHO<Dim>* lowerPb = dynamic_cast<P_LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		FunctionalBasis<Dim - 1>* higherBasis = higherPb->HHO->FaceBasis;
		FunctionalBasis<Dim - 1>* lowerBasis = lowerPb->HHO->FaceBasis;
		int nHigherUnknowns = higherBasis->Size();
		int nLowerUnknowns = lowerBasis->Size();

		assert(higherBasis->BasisCode().compare(lowerBasis->BasisCode()) == 0 && "The bases must be the same");

		if (!higherBasis->IsHierarchical)
			Utils::FatalError("The natural injection is not implemented for non-hierarchical bases.");

		FaceParallelLoop<Dim> parallelLoop(mesh->Faces);
		parallelLoop.ReserveChunkCoeffsSize(lowerBasis->LocalFunctions.size());
		parallelLoop.Execute([this, nHigherUnknowns, nLowerUnknowns](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				if (f->HasDirichletBC())
					return;

				for (int i = 0; i < nLowerUnknowns; i++)
					chunk->Results.Coeffs.Add(f->Number*nHigherUnknowns + i, f->Number*nLowerUnknowns + i, 1);
			});

		P = SparseMatrix(higherPb->HHO->nInteriorAndNeumannFaces * nHigherUnknowns, lowerPb->HHO->nInteriorAndNeumannFaces * nLowerUnknowns);
		parallelLoop.Fill(P);
	}

public:
	void SetupRestriction() override
	{
		R = P.transpose().eval();
	}
};



template <int Dim>
class P_MultigridForHHO : public Multigrid
{
private:
	Diffusion_HHO<Dim>* _problem;
public:

	P_MultigridForHHO(Diffusion_HHO<Dim>* problem)
		: Multigrid(0)
	{
		this->HP_CS = HP_CoarsStgy::P_only;
		this->_problem = problem;
		this->BlockSizeForBlockSmoothers = problem->HHO->nFaceUnknowns;
	}

	void BeginSerialize(ostream& os) const override
	{
		os << "p-MultigridForHHO" << endl;
		os << "\t" << "Prolongation            : natural injection";
		os << endl;
	}

protected:
	Level* CreateFineLevel() const override
	{
		return new P_LevelForHHO<Dim>(0, _problem);
	}

	Level* CreateCoarseLevel(Level* fineLevel, CoarseningType coarseningType) override
	{
		if (coarseningType != CoarseningType::P)
			Utils::FatalError("Only p-coarsening allowed for this multigrid.");

		P_LevelForHHO<Dim>* hhoFineLevel = dynamic_cast<P_LevelForHHO<Dim>*>(fineLevel);
		Diffusion_HHO<Dim>* lowerDegreeProblem = hhoFineLevel->_problem->GetProblemForLowerDegree();
		P_LevelForHHO<Dim>* coarseLevel = new P_LevelForHHO<Dim>(fineLevel->Number + 1, lowerDegreeProblem);
		return coarseLevel;
	}

	void InitializeCoarseSolver() override
	{
		MultigridForHHO<Dim>* mg = dynamic_cast<MultigridForHHO<Dim>*>(this->CoarseSolver);
		if (mg)
		{
			Level* level = this->_fineLevel;
			while (level->CoarserLevel != nullptr)
				level = level->CoarserLevel;

			P_LevelForHHO<Dim>* lowOrderLevel = dynamic_cast<P_LevelForHHO<Dim>*>(level);
			Diffusion_HHO<Dim>* lowOrderPb = lowOrderLevel->_problem;

			mg->InitializeWithProblem(lowOrderPb);
		}
	}
};