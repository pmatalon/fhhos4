#pragma once
#include <Eigen/Sparse>
#include <vector>
#include "BlockSOR.h"
#include "Multigrid.h"
#include "../Mesh/Mesh.h"
#include "../HHO/Poisson_HHO_Element.h"
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

		auto I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
		auto J_f_c = GetGlobalCanonicalInjectionMatrixCoarseToFine();
		//auto Pi_c = GetGlobalProjectorMatrixFromCellsToFaces(coarsePb); // to be removed later
		auto Pi_f = GetGlobalProjectorMatrixFromCellsToFaces(finePb);

		/*finePb->ExportMatrix(I_c, "I_c");
		finePb->ExportMatrix(J_f_c, "J_f_c");
		finePb->ExportMatrix(Pi_c, "Pi_c"); // to be removed later
		finePb->ExportMatrix(Pi_f, "Pi_f");*/

		P = (Pi_f * J_f_c * I_c).pruned();

		finePb->ExportMatrix(P, "P");
	}

	void SetupRestriction()
	{
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

		Eigen::SparseMatrix<double> M_c = GetGlobalMassMatrix(coarsePb); // Coarse mass matrix
		Eigen::SparseMatrix<double> M_f = GetGlobalMassMatrix(finePb); // Fine mass matrix

		finePb->ExportMatrix(M_c, "level_" + to_string(this->Number) + "_M_c_qui_marche_pas");

		Eigen::SparseLU<Eigen::SparseMatrix<double>> solverM_c;
		solverM_c.analyzePattern(M_c);
		solverM_c.factorize(M_c);

		//R = invM_c * P.transpose() * M_f;
		R = solverM_c.solve(P.transpose() * M_f);

		finePb->ExportMatrix(R, "R");
	}

	Eigen::SparseMatrix<double> GetGlobalMassMatrix(Poisson_HHO<Dim>* problem)
	{
		int faceLocalUnknowns = problem->HHO.nLocalFaceUnknowns;

		ParallelLoop<Element<Dim>*> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(4 * faceLocalUnknowns * faceLocalUnknowns);
		parallelLoop.Execute([this, faceLocalUnknowns](Element<Dim>* e, ParallelChunk* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);
				for (auto f : element->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

					double weight = element->Measure() / element->FrontierMeasure();
					Eigen::MatrixXd M_face = face->FaceMassMatrix();

					BigNumber faceGlobalNumber = face->Number;
					BigNumber faceLocalNumber = element->LocalNumberOf(face);

					chunk->Coeffs.Add(faceGlobalNumber*faceLocalUnknowns, faceGlobalNumber*faceLocalUnknowns, weight * M_face);
				}
			});
		Eigen::SparseMatrix<double> M(problem->HHO.nTotalFaceUnknowns, problem->HHO.nTotalFaceUnknowns);
		parallelLoop.Fill(M);

		return M;
	}

	Eigen::SparseMatrix<double> GetGlobalInterpolationMatrixFromFacesToCells(Poisson_HHO<Dim>* problem)
	{
		int nCellUnknowns = problem->HHO.nLocalCellUnknowns;
		int nFaceUnknowns = problem->HHO.nLocalFaceUnknowns;

		ParallelLoop<Element<Dim>*> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nCellUnknowns, nFaceUnknowns](Element<Dim>* e, ParallelChunk* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);

				Eigen::MatrixXd M_cell = element->CellMassMatrix();
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
					Eigen::MatrixXd Pi_face = face->GetProjFromCell(element);

					double weight = element->Measure() / element->FrontierMeasure();
					chunk->Coeffs.Add(elemGlobalNumber*nCellUnknowns, faceGlobalNumber*nFaceUnknowns, weight * invM_cell * Pi_face.transpose() * M_face);
				}
			});
		Eigen::SparseMatrix<double> I(problem->HHO.nTotalCellUnknowns, problem->HHO.nTotalFaceUnknowns);
		parallelLoop.Fill(I);

		return I;
	}

	Eigen::SparseMatrix<double> GetGlobalProjectorMatrixFromCellsToFaces(Poisson_HHO<Dim>* problem)
	{
		int nCellUnknowns = problem->HHO.nLocalCellUnknowns;
		int nFaceUnknowns = problem->HHO.nLocalFaceUnknowns;

		ParallelLoop<Element<Dim>*> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nCellUnknowns, nFaceUnknowns](Element<Dim>* e, ParallelChunk* chunk)
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

					double weight = element->Measure() / element->FrontierMeasure();
					chunk->Coeffs.Add(faceGlobalNumber*nFaceUnknowns, elemGlobalNumber*nCellUnknowns, weight*face->GetProjFromCell(element));
				}
			});

		Eigen::SparseMatrix<double> Pi(problem->HHO.nTotalFaceUnknowns, problem->HHO.nTotalCellUnknowns);
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

		int nCoarseFaceUnknowns = coarsePb->HHO.nLocalFaceUnknowns;
		int nCoarseCellUnknowns = coarsePb->HHO.nLocalCellUnknowns;

		int nFineFaceUnknowns = finePb->HHO.nLocalFaceUnknowns;
		int nFineCellUnknowns = finePb->HHO.nLocalCellUnknowns;

		ParallelLoop<Element<Dim>*> parallelLoop(coarseMesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nFineCellUnknowns * 4 * nCoarseCellUnknowns);

		parallelLoop.Execute([this, nCoarseFaceUnknowns, nCoarseCellUnknowns, nFineFaceUnknowns, nFineCellUnknowns](Element<Dim>* ce, ParallelChunk* chunk)
			{
				Poisson_HHO_Element<Dim>* coarseElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(ce);

				Eigen::MatrixXd local_J_f_c = coarseElement->ComputeCanonicalInjectionMatrixCoarseToFine();
				for (auto fineElement : coarseElement->FinerElements)
				{
					BigNumber coarseElemGlobalNumber = coarseElement->Number;
					BigNumber fineElemGlobalNumber = fineElement->Number;
					BigNumber fineElemLocalNumber = coarseElement->LocalNumberOf(fineElement);

					chunk->Coeffs.Add(fineElemGlobalNumber*nFineCellUnknowns, coarseElemGlobalNumber*nCoarseCellUnknowns, local_J_f_c.block(fineElemLocalNumber*nFineCellUnknowns, 0, nFineCellUnknowns, nCoarseCellUnknowns));
				}
			});

		Eigen::SparseMatrix<double> J_f_c(finePb->HHO.nTotalCellUnknowns, coarsePb->HHO.nTotalCellUnknowns);
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
	//const int DefaultMatrixMaxSizeForCoarsestLevel = 100;
public:
	int MatrixMaxSizeForCoarsestLevel;

	MultigridForHHO(Poisson_HHO<Dim>* problem) : MultigridForHHO(problem, 0)
	{}
	//MultigridForHHO(Poisson_HHO<Dim>* problem, int nLevels) : MultigridForHHO(problem, nLevels <= 0, nLevels > 0 ? nLevels : DefaultMatrixMaxSizeForCoarsestLevel)
	//{}

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