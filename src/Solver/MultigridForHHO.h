#pragma once
#include <Eigen/Sparse>
#include <vector>
#include "BlockSOR.h"
#include "Multigrid.h"
#include "../Mesh/Mesh.h"
#include "../HHO/Poisson_HHO_Element.h"
#include "../Utils/ElementParallelLoop.h"
using namespace std;

template <int Dim>
class LevelForHHO : public Level
{
private:
	Poisson_HHO<Dim>* _problem;
public:
	LevelForHHO(int number, Poisson_HHO<Dim>* problem, bool useGalerkinOperator)
		: Level(number, new BlockGaussSeidelSmoother(problem->HHO.nLocalFaceUnknowns, 1), new ReverseBlockGaussSeidelSmoother(problem->HHO.nLocalFaceUnknowns, 1))
	{
		this->UseGalerkinOperator = useGalerkinOperator;
		this->_problem = problem;
	}

	void ExportVector(Eigen::VectorXd& v, string suffix)
	{
		this->_problem->ExportVector(v, suffix);
	}

private:
	void SetupDiscretizedOperator() override 
	{
		this->OperatorMatrix = this->_problem->A;
	}

	void OnStartSetup() override
	{
		cout << "\t\tMesh              : " << this->_problem->_mesh->Elements.size() << " elements" << endl;
		if (!IsCoarsestLevel())
		{
			Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

			if (this->UseGalerkinOperator)
				coarsePb->InitHHO();
			else
				coarsePb->Assemble(Action::None);
		}
	}

	void SetupProlongation() override
	{
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

		SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
		SparseMatrix J_f_c = GetGlobalCanonicalInjectionMatrixCoarseToFine();
		SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCellsToFaces(finePb);
		
		finePb->ExportMatrix(I_c, "I_c");
		finePb->ExportMatrix(J_f_c, "J_f_c");
		finePb->ExportMatrix(Pi_f, "Pi_f");

		P = Pi_f * J_f_c * I_c;

		finePb->ExportMatrix(P, "P");
	}

	void SetupRestriction() override
	{
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

		SparseMatrix invM_c = GetInverseGlobalMassMatrix_Faces(coarsePb); // Coarse mass matrix
		SparseMatrix M_f = GetGlobalMassMatrix_Faces(finePb); // Fine mass matrix

		//finePb->ExportMatrix(invM_c, "invM_c");
		//finePb->ExportMatrix(M_f, "M_f");

		R = invM_c * P.transpose() * M_f;

		finePb->ExportMatrix(R, "R");
	}

	inline double Weight(Element<Dim>* element, Face<Dim>* face)
	{
		return element->Measure() / (face->Element1->Measure() + face->Element2->Measure());
	}

	SparseMatrix GetGlobalMassMatrix_Faces(Poisson_HHO<Dim>* problem)
	{
		int faceLocalUnknowns = problem->HHO.nLocalFaceUnknowns;

		FaceParallelLoop<Dim> parallelLoop(problem->_mesh->InteriorFaces);
		parallelLoop.ReserveChunkCoeffsSize(faceLocalUnknowns * faceLocalUnknowns);
		parallelLoop.Execute([this, faceLocalUnknowns](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
				chunk->Results.Coeffs.Add(face->Number*faceLocalUnknowns, face->Number*faceLocalUnknowns, face->FaceMassMatrix() / face->Measure());
			});
		SparseMatrix M(problem->HHO.nTotalFaceUnknowns, problem->HHO.nTotalFaceUnknowns);
		parallelLoop.Fill(M);

		return M;
	}

	SparseMatrix GetInverseGlobalMassMatrix_Faces(Poisson_HHO<Dim>* problem)
	{
		int faceLocalUnknowns = problem->HHO.nLocalFaceUnknowns;

		FaceParallelLoop<Dim> parallelLoop(problem->_mesh->InteriorFaces);
		parallelLoop.ReserveChunkCoeffsSize(faceLocalUnknowns * faceLocalUnknowns);
		parallelLoop.Execute([this, faceLocalUnknowns](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
				chunk->Results.Coeffs.Add(face->Number*faceLocalUnknowns, face->Number*faceLocalUnknowns, face->InvFaceMassMatrix() * face->Measure());
			});
		SparseMatrix M(problem->HHO.nTotalFaceUnknowns, problem->HHO.nTotalFaceUnknowns);
		parallelLoop.Fill(M);

		return M;
	}

	SparseMatrix GetInverseGlobalMassMatrix_Cells(Poisson_HHO<Dim>* problem)
	{
		int cellLocalUnknowns = problem->HHO.nLocalCellUnknowns;

		ParallelLoop<Element<Dim>*> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(cellLocalUnknowns * cellLocalUnknowns);
		parallelLoop.Execute([this, cellLocalUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);
				Eigen::MatrixXd M_cell = element->CellMassMatrix();
				Eigen::MatrixXd invM_cell = M_cell.inverse();
				chunk->Results.Coeffs.Add(element->Number*cellLocalUnknowns, element->Number*cellLocalUnknowns, invM_cell);
			});
		SparseMatrix M(problem->HHO.nTotalCellUnknowns, problem->HHO.nTotalCellUnknowns);
		parallelLoop.Fill(M);

		return M;
	}

	FunctionalBasis<Dim>* GetCellInterpolationBasis(Poisson_HHO<Dim>* problem)
	{
		FunctionalBasis<Dim-1>* faceBasis = dynamic_cast<Poisson_HHO_Element<Dim>*>(problem->_mesh->Elements[0])->FaceBasis;
		int faceDegreePSpace = faceBasis->GetDegree();

		//int cellDegreeQSpace = (int)floor(2 * sqrt(faceDegreePSpace + 1) - 1);
		//FunctionalBasis<Dim>* cellInterpolationBasis = new FunctionalBasis<Dim>("monomials", cellDegreeQSpace, true);
		int cellDegreePSpace = faceDegreePSpace + 1;
		FunctionalBasis<Dim>* cellInterpolationBasis = new FunctionalBasis<Dim>(faceBasis->BasisCode(), cellDegreePSpace, false);

		return cellInterpolationBasis;
	}

	SparseMatrix GetGlobalInterpolationMatrixFromFacesToCells(Poisson_HHO<Dim>* problem)
	{
		//return GetAdjointToGlobalProjUsingCellInnerProduct(problem);
		return GetReconstructionMatrix(problem);
	}



	/*SparseMatrix GetAdjointToGlobalProjUsingCellInnerProduct(Poisson_HHO<Dim>* problem)
	{
		int nCellUnknowns = problem->HHO.nLocalCellUnknowns;
		int nFaceUnknowns = problem->HHO.nLocalFaceUnknowns;

		FunctionalBasis<Dim>* cellInterpolationBasis = GetCellInterpolationBasis(problem);

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nCellUnknowns, nFaceUnknowns, cellInterpolationBasis](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);

				Eigen::MatrixXd M_cell = element->CellMassMatrix(cellInterpolationBasis);
				Eigen::MatrixXd invM_cell = M_cell.inverse();

				for (auto f : element->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

					BigNumber elemGlobalNumber = element->Number;
					BigNumber faceGlobalNumber = face->Number;
					BigNumber faceLocalNumber = element->LocalNumberOf(face);

					Eigen::MatrixXd M_face = face->FaceMassMatrix();
					Eigen::MatrixXd Pi_face = face->GetProjFromCell(element, cellInterpolationBasis);

					double weight = Weight(element, face);
					chunk->Results.Coeffs.Add(elemGlobalNumber*cellInterpolationBasis->Size(), faceGlobalNumber*nFaceUnknowns, weight * invM_cell * Pi_face.transpose() * M_face);
				}
			});

		SparseMatrix I(problem->HHO.nElements * cellInterpolationBasis->Size(), problem->HHO.nTotalFaceUnknowns);
		parallelLoop.Fill(I);
		delete cellInterpolationBasis;

		return I;
	}*/

	SparseMatrix GetReconstructionMatrix(Poisson_HHO<Dim>* problem)
	{
		int nCellUnknowns = problem->HHO.nLocalCellUnknowns;
		int nFaceUnknowns = problem->HHO.nLocalFaceUnknowns;

		FunctionalBasis<Dim>* cellInterpolationBasis = GetCellInterpolationBasis(problem);

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nFaceUnknowns, cellInterpolationBasis](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);
				for (auto f : element->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

					BigNumber elemGlobalNumber = element->Number;
					BigNumber faceGlobalNumber = face->Number;
					BigNumber faceLocalNumber = element->LocalNumberOf(face);
					Eigen::MatrixXd reconstructMatrix = element->ReconstructionFromFacesMatrix();
					chunk->Results.Coeffs.Add(elemGlobalNumber * cellInterpolationBasis->Size(), faceGlobalNumber * nFaceUnknowns, reconstructMatrix.block(0, faceLocalNumber*nFaceUnknowns, cellInterpolationBasis->Size(), nFaceUnknowns));
				}
			});

		SparseMatrix M(problem->HHO.nElements * cellInterpolationBasis->Size(), problem->HHO.nTotalFaceUnknowns);
		parallelLoop.Fill(M);
		delete cellInterpolationBasis;

		return M;
	}

	SparseMatrix GetGlobalProjectorMatrixFromCellsToFaces(Poisson_HHO<Dim>* problem)
	{
		int nFaceUnknowns = problem->HHO.nLocalFaceUnknowns;

		FunctionalBasis<Dim>* cellInterpolationBasis = GetCellInterpolationBasis(problem);
		int nCellUnknowns = cellInterpolationBasis->Size();

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nCellUnknowns, nFaceUnknowns, cellInterpolationBasis](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);

				for (auto f : element->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

					BigNumber elemGlobalNumber = element->Number;
					BigNumber faceGlobalNumber = face->Number;
					BigNumber faceLocalNumber = element->LocalNumberOf(face);

					double weight = Weight(element, face);
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, elemGlobalNumber*nCellUnknowns, weight*face->GetProjFromCell(element, cellInterpolationBasis));
				}
			});

		SparseMatrix Pi(problem->HHO.nTotalFaceUnknowns, problem->HHO.nElements * nCellUnknowns);
		parallelLoop.Fill(Pi);

		return Pi;
	}

	SparseMatrix GetGlobalCanonicalInjectionMatrixCoarseToFine()
	{
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		Mesh<Dim>* coarseMesh = coarsePb->_mesh;

		FunctionalBasis<Dim>* cellInterpolationBasis = GetCellInterpolationBasis(finePb);
		int nCellUnknowns = cellInterpolationBasis->Size();

		ElementParallelLoop<Dim> parallelLoop(coarseMesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nCellUnknowns);

		parallelLoop.Execute([this, cellInterpolationBasis, nCellUnknowns](Element<Dim>* ce, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* coarseElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(ce);

				Eigen::MatrixXd local_J_f_c = coarseElement->ComputeCanonicalInjectionMatrixCoarseToFine(cellInterpolationBasis);
				for (auto fineElement : coarseElement->FinerElements)
				{
					BigNumber coarseElemGlobalNumber = coarseElement->Number;
					BigNumber fineElemGlobalNumber = fineElement->Number;
					BigNumber fineElemLocalNumber = coarseElement->LocalNumberOf(fineElement);

					chunk->Results.Coeffs.Add(fineElemGlobalNumber*nCellUnknowns, coarseElemGlobalNumber*nCellUnknowns, local_J_f_c.block(fineElemLocalNumber*nCellUnknowns, 0, nCellUnknowns, nCellUnknowns));
				}
			});

		SparseMatrix J_f_c(finePb->HHO.nElements * nCellUnknowns, coarsePb->HHO.nElements * nCellUnknowns);
		parallelLoop.Fill(J_f_c);
		return J_f_c;
	}

public:
	~LevelForHHO()
	{
		if (!this->IsFinestLevel())
			delete _problem;
	}
};

template <int Dim>
class MultigridForHHO : public Multigrid
{
private:
	Poisson_HHO<Dim>* _problem;
	bool _automaticNumberOfLevels;
	int _nLevels;
public:
	int MatrixMaxSizeForCoarsestLevel;
	bool UseGalerkinOperator;

	MultigridForHHO(Poisson_HHO<Dim>* problem) : MultigridForHHO(problem, 0)
	{}

	MultigridForHHO(Poisson_HHO<Dim>* problem, int nLevels) : Multigrid()
	{
		this->_problem = problem;
		this->_automaticNumberOfLevels = nLevels <= 0;
		
		if (this->_automaticNumberOfLevels)
		{
			this->_nLevels = 0;
			this->MatrixMaxSizeForCoarsestLevel = 100;
		}
		else
		{
			this->_nLevels = nLevels;
			this->MatrixMaxSizeForCoarsestLevel = 0;
		}

		this->_fineLevel = new LevelForHHO<Dim>(0, problem, false);
	}
	
	void Serialize(ostream& os) const override
	{
		os << "MultigridForHHO" << endl;
		os << "\t" << "Levels        : ";
		if (_automaticNumberOfLevels && _nLevels == 0)
			os << "automatic coarsening until matrix size <= " << MatrixMaxSizeForCoarsestLevel << endl;
		else if (_automaticNumberOfLevels && _nLevels > 0)
			os << _nLevels << " (automatic)" << endl;
		else
			os << _nLevels << endl;
		os << "\t" << "Pre-smoothing : " << *(_fineLevel->PreSmoother) << endl;
		os << "\t" << "Post-smoothing: " << *(_fineLevel->PostSmoother);
	}

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);

		cout << "Setup..." << endl;

		CoarseningStrategy coarseningStrategy = CoarseningStrategy::AgglomerationAndMergeColinearFaces;
		//CoarseningStrategy coarseningStrategy = CoarseningStrategy::AgglomerationAndKeepFineFaces;

		LevelForHHO<Dim>* finerLevel = dynamic_cast<LevelForHHO<Dim>*>(this->_fineLevel);
		finerLevel->OperatorMatrix = A;
		Poisson_HHO<Dim>* problem = _problem;
		if (!_automaticNumberOfLevels)
		{
			for (int levelNumber = 1; levelNumber < _nLevels; levelNumber++)
			{
				//cout << "fine mesh" << endl << *(problem->_mesh) << endl << endl;
				problem->_mesh->CoarsenMesh(coarseningStrategy);
				if (problem->_mesh->CoarseMesh->InteriorFaces.size() == 0)
				{
					cout << "Warning: impossible to coarsen the mesh any more.";
					break;
				}
				problem = problem->GetProblemOnCoarserMesh();
				//cout << "coarse mesh" << endl << *(problem->_mesh) << endl << endl;
				LevelForHHO<Dim>* coarseLevel = new LevelForHHO<Dim>(levelNumber, problem, UseGalerkinOperator);

				finerLevel->CoarserLevel = coarseLevel;
				coarseLevel->FinerLevel = finerLevel;

				finerLevel->Setup();

				finerLevel = coarseLevel;
			}
			finerLevel->Setup();
			this->SetupCoarseSolver();
		}
		else
		{
			int levelNumber = 0;
			while (problem->HHO.nTotalFaceUnknowns > MatrixMaxSizeForCoarsestLevel)
			{
				problem->_mesh->CoarsenMesh(coarseningStrategy);
				auto nCoarseInteriorFaces = problem->_mesh->CoarseMesh->InteriorFaces.size();
				auto nCoarseUnknowns = nCoarseInteriorFaces * problem->HHO.nLocalFaceUnknowns;
				if (nCoarseInteriorFaces == 0)
					break;
				problem = problem->GetProblemOnCoarserMesh();
				levelNumber++;
				LevelForHHO<Dim>* coarseLevel = new LevelForHHO<Dim>(levelNumber, problem, UseGalerkinOperator);

				finerLevel->CoarserLevel = coarseLevel;
				coarseLevel->FinerLevel = finerLevel;

				finerLevel->Setup();

				finerLevel = coarseLevel;
			}
			finerLevel->Setup();
			_nLevels = levelNumber + 1;
			this->SetupCoarseSolver();
			cout << "\t--> " << _nLevels << " levels built." << endl;
		}
	}

	//~MultigridForHHO()
	//{}
};