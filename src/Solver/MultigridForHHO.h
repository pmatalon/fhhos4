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
	Eigen::SparseMatrix<double> R;
	Eigen::SparseMatrix<double> P;
	Poisson_HHO<Dim>* _problem;
public:
	LevelForHHO(int number, Poisson_HHO<Dim>* problem)
		: Level(number, new BlockGaussSeidelSmoother(problem->HHO.nLocalFaceUnknowns, 1), new ReverseBlockGaussSeidelSmoother(problem->HHO.nLocalFaceUnknowns, 1))
	{
		this->_problem = problem;
	}

	void Setup()
	{
		cout << "\tSetup level " << this->Number << endl;
		/*this->OperatorMatrix = this->_problem->A;
		if (!this->IsCoarsestLevel())
		{
			SetupProlongation();
			SetupRestriction();
			SetupSmoothers();
		}*/

		// Galerkin operator
		if (!this->IsCoarsestLevel())
		{
			cout << "\t\tProlongation      : "; cout.flush();
			SetupProlongation();
			cout << Utils::MatrixInfo(this->P, "P") << endl;

			cout << "\t\tRestriction       : "; cout.flush();
			SetupRestriction();
			cout << Utils::MatrixInfo(this->R, "R") << endl;
		}
		
		if (this->IsFinestLevel())
		{
			cout << "\t\tFine grid operator: "; cout.flush();
			this->OperatorMatrix = this->_problem->A;
			cout << Utils::MatrixInfo(this->OperatorMatrix, "A") << endl;
		}
		else
		{
			LevelForHHO<Dim>* finerLevel = dynamic_cast<LevelForHHO<Dim>*>(FinerLevel);
			cout << "\t\tGalerkin operator : "; cout.flush();
			this->OperatorMatrix = (finerLevel->R * finerLevel->OperatorMatrix * finerLevel->P).pruned();
			cout << Utils::MatrixInfo(this->OperatorMatrix, "A") << endl;
		}

		if (!this->IsCoarsestLevel())
		{
			cout << "\t\tSmoothers..." << endl;
			SetupSmoothers();
		}		
	}

	Eigen::VectorXd Restrict(Eigen::VectorXd& vectorOnThisLevel) override
	{
		Eigen::VectorXd coarseVector = R * vectorOnThisLevel;
		return coarseVector;
	}
	
	Eigen::VectorXd Prolong(Eigen::VectorXd& vectorOnTheCoarserLevel)  override
	{
		Eigen::VectorXd vectorOnThisLevel = P * vectorOnTheCoarserLevel;
		return vectorOnThisLevel;
	}

	void ExportVector(Eigen::VectorXd& v, string suffix)
	{
		this->_problem->ExportVector(v, suffix);
	}

private:

	void SetupProlongation()
	{
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

		//Eigen::SparseMatrix<double> Pi_c = GetGlobalProjectorMatrixFromCellsToFaces(coarsePb); // to be removed later
		Eigen::SparseMatrix<double> I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
		/*Eigen::SparseMatrix<double> invM_cells = GetInverseGlobalMassMatrix_Cells(coarsePb);
		Eigen::SparseMatrix<double> M_faces = GetGlobalMassMatrix_Faces(coarsePb);
		Eigen::SparseMatrix<double> I_c = invM_cells * Pi_c.transpose() * M_faces;*/
		//Eigen::SparseMatrix<double> M_faces_f = GetGlobalMassMatrix_Faces(finePb); // to be removed later
		//finePb->ExportMatrix(M_faces_f, "M_faces_f"); // to be removed later

		Eigen::SparseMatrix<double> J_f_c = GetGlobalCanonicalInjectionMatrixCoarseToFine();
		Eigen::SparseMatrix<double> Pi_f = GetGlobalProjectorMatrixFromCellsToFaces(finePb);
		
		//cout << "----------------- Pi_f (" << Pi_f.rows() << ", " << Pi_f.cols() << ") ----------------" << endl << Pi_f << endl;
		finePb->ExportMatrix(I_c, "I_c");
		finePb->ExportMatrix(J_f_c, "J_f_c");
		//finePb->ExportMatrix(Pi_c, "Pi_c"); // to be removed later
		finePb->ExportMatrix(Pi_f, "Pi_f");

		P = (Pi_f * J_f_c * I_c).pruned();

		finePb->ExportMatrix(P, "P");
	}

	void SetupRestriction()
	{
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

		Eigen::SparseMatrix<double> invM_c = GetInverseGlobalMassMatrix_Faces(coarsePb); // Coarse mass matrix
		Eigen::SparseMatrix<double> M_f = GetGlobalMassMatrix_Faces(finePb); // Fine mass matrix

		//finePb->ExportMatrix(M_f, "level_" + to_string(this->Number) + "_M_f");

		/*Eigen::SparseLU<Eigen::SparseMatrix<double>> solverM_c;
		solverM_c.analyzePattern(M_c);
		solverM_c.factorize(M_c);*/

		//R = solverM_c.solve(P.transpose() * M_f);
		R = invM_c * P.transpose() * M_f;

		finePb->ExportMatrix(R, "R");
	}

	inline double Weight(Element<Dim>* element, Face<Dim>* face)
	{
		return 0.5;//element->Measure() / element->FrontierMeasure();
	}

	Eigen::SparseMatrix<double> GetGlobalMassMatrix_Faces(Poisson_HHO<Dim>* problem)
	{
		int faceLocalUnknowns = problem->HHO.nLocalFaceUnknowns;

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(4 * faceLocalUnknowns * faceLocalUnknowns);
		parallelLoop.Execute([this, faceLocalUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);
				for (auto f : element->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

					double weight = Weight(element, face);
					Eigen::MatrixXd M_face = face->FaceMassMatrix();

					BigNumber faceGlobalNumber = face->Number;
					BigNumber faceLocalNumber = element->LocalNumberOf(face);

					chunk->Results.Coeffs.Add(faceGlobalNumber*faceLocalUnknowns, faceGlobalNumber*faceLocalUnknowns, weight * M_face);
				}
			});
		Eigen::SparseMatrix<double> M(problem->HHO.nTotalFaceUnknowns, problem->HHO.nTotalFaceUnknowns);
		parallelLoop.Fill(M);

		return M;
	}

	Eigen::SparseMatrix<double> GetInverseGlobalMassMatrix_Faces(Poisson_HHO<Dim>* problem)
	{
		int faceLocalUnknowns = problem->HHO.nLocalFaceUnknowns;

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(4 * faceLocalUnknowns * faceLocalUnknowns);
		parallelLoop.Execute([this, faceLocalUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);
				for (auto f : element->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

					double weight = Weight(element, face);
					Eigen::MatrixXd invM_face = face->InvFaceMassMatrix();

					BigNumber faceGlobalNumber = face->Number;
					BigNumber faceLocalNumber = element->LocalNumberOf(face);

					chunk->Results.Coeffs.Add(faceGlobalNumber*faceLocalUnknowns, faceGlobalNumber*faceLocalUnknowns, invM_face / weight);
				}
			});
		Eigen::SparseMatrix<double> M(problem->HHO.nTotalFaceUnknowns, problem->HHO.nTotalFaceUnknowns);
		parallelLoop.Fill(M);

		return M;
	}

	Eigen::SparseMatrix<double> GetInverseGlobalMassMatrix_Cells(Poisson_HHO<Dim>* problem)
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
		Eigen::SparseMatrix<double> M(problem->HHO.nTotalCellUnknowns, problem->HHO.nTotalCellUnknowns);
		parallelLoop.Fill(M);

		return M;
	}

	FunctionalBasis<Dim>* GetCellInterpolationBasis(Poisson_HHO<Dim>* problem)
	{
		int faceDegreePSpace = dynamic_cast<Poisson_HHO_Element<Dim>*>(problem->_mesh->Elements[0])->FaceBasis->GetDegree();
		int cellDegreeQSpace = (int)floor(2 * sqrt(faceDegreePSpace + 1) - 1);
		FunctionalBasis<Dim>* cellInterpolationBasis = new FunctionalBasis<Dim>("monomials", cellDegreeQSpace, true);
		return cellInterpolationBasis;
	}

	Eigen::SparseMatrix<double> GetGlobalInterpolationMatrixFromFacesToCells(Poisson_HHO<Dim>* problem)
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

		Eigen::SparseMatrix<double> I(problem->HHO.nElements * cellInterpolationBasis->Size(), problem->HHO.nTotalFaceUnknowns);
		parallelLoop.Fill(I);
		delete cellInterpolationBasis;

		return I;
	}

	Eigen::SparseMatrix<double> GetGlobalProjectorMatrixFromCellsToFaces(Poisson_HHO<Dim>* problem)
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

		Eigen::SparseMatrix<double> Pi(problem->HHO.nTotalFaceUnknowns, problem->HHO.nElements * nCellUnknowns);
		parallelLoop.Fill(Pi);

		/*NonZeroCoefficients Pi_coeffs(problem->HHO.nTotalCellUnknowns * problem->HHO.nTotalFaceUnknowns);
		int nCellUnknowns = problem->HHO.nLocalCellUnknowns;
		int nFaceUnknowns = problem->HHO.nLocalFaceUnknowns;
		for (auto f : problem->_mesh->Faces)
		{
			if (f->IsDomainBoundary)
				continue;

			Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

			Poisson_HHO_Element<Dim>* element1 = dynamic_cast<Poisson_HHO_Element<Dim>*>(face->Element1);
			Poisson_HHO_Element<Dim>* element2 = dynamic_cast<Poisson_HHO_Element<Dim>*>(face->Element2);
			//Eigen::MatrixXd local_Pi = element->ComputeProjectorMatrixFromCellOntoFaces();

			//cout << "Face " << face->Number << " takes value from element " << element->Number << endl;

			BigNumber faceGlobalNumber = face->Number;
			BigNumber elem1GlobalNumber = element1->Number;
			BigNumber faceLocalNumber1 = element1->LocalNumberOf(face);
			BigNumber elem2GlobalNumber = element2->Number;
			BigNumber faceLocalNumber2 = element2->LocalNumberOf(face);

			Pi_coeffs.Add(faceGlobalNumber*nFaceUnknowns, elem1GlobalNumber*nCellUnknowns, 0.5 * face->GetProjFromCell(element1));//local_Pi.block(faceLocalNumber*nFaceUnknowns, 0, nFaceUnknowns, nCellUnknowns));
			Pi_coeffs.Add(faceGlobalNumber*nFaceUnknowns, elem2GlobalNumber*nCellUnknowns, 0.5 * face->GetProjFromCell(element2));
		}

		Eigen::SparseMatrix<double> Pi(problem->HHO.nTotalFaceUnknowns, problem->HHO.nTotalCellUnknowns);
		Pi_coeffs.Fill(Pi);*/
		return Pi;
	}

	Eigen::SparseMatrix<double> GetGlobalCanonicalInjectionMatrixCoarseToFine()
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

		Eigen::SparseMatrix<double> J_f_c(finePb->HHO.nElements * nCellUnknowns, coarsePb->HHO.nElements * nCellUnknowns);
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

		this->_fineLevel = new LevelForHHO<Dim>(0, problem);
	}
	
	virtual void Serialize(ostream& os) const override
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

	void Setup(const Eigen::SparseMatrix<double>& A) override
	{
		IterativeSolver::Setup(A);

		cout << "Setup..." << endl;

		LevelForHHO<Dim>* finerLevel = dynamic_cast<LevelForHHO<Dim>*>(this->_fineLevel);
		Poisson_HHO<Dim>* problem = _problem;
		if (!_automaticNumberOfLevels)
		{
			for (int levelNumber = 1; levelNumber < _nLevels; levelNumber++)
			{
				dynamic_cast<CartesianGrid2D*>(problem->_mesh)->BuildCoarserMesh();
				problem = problem->GetProblemOnCoarserMesh();
				problem->Assemble(Action::None);
				LevelForHHO<Dim>* coarseLevel = new LevelForHHO<Dim>(levelNumber, problem);

				finerLevel->CoarserLevel = coarseLevel;
				coarseLevel->FinerLevel = finerLevel;

				finerLevel->Setup();

				finerLevel = coarseLevel;
			}
			finerLevel->Setup();
		}
		else
		{
			int levelNumber = 0;
			while (problem->A.rows() > MatrixMaxSizeForCoarsestLevel)
			{
				levelNumber++;
				dynamic_cast<CartesianGrid2D*>(problem->_mesh)->BuildCoarserMesh();
				problem = problem->GetProblemOnCoarserMesh();
				problem->Assemble(Action::None);
				LevelForHHO<Dim>* coarseLevel = new LevelForHHO<Dim>(levelNumber, problem);

				finerLevel->CoarserLevel = coarseLevel;
				coarseLevel->FinerLevel = finerLevel;

				finerLevel->Setup();

				finerLevel = coarseLevel;
			}
			finerLevel->Setup();
			_nLevels = levelNumber + 1;
			cout << "\t--> " << _nLevels << " levels built." << endl;
		}
	}

	//~MultigridForHHO()
	//{}
};