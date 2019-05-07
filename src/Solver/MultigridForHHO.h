#pragma once
#include <Eigen/Sparse>
#include <vector>
#include "GeometricMultigrid.h"
#include "BlockSOR.h"
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
		: Level(number, Smoother(new BlockGaussSeidel(problem->HHO.nLocalFaceUnknowns), 1), Smoother(new ReverseBlockGaussSeidel(problem->HHO.nLocalFaceUnknowns), 1))
	{
		this->_problem = problem;
	}

	void Setup()
	{
		this->OperatorMatrix = this->_problem->A;
		if (!this->IsCoarsestLevel())
		{
			SetupProlongation();
			SetupRestriction();
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
		//P = Eigen::SparseMatrix<double>(this->OperatorMatrix.rows(), this->CoarserLevel->OperatorMatrix.rows());

		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		Mesh<Dim>* coarseMesh = coarsePb->_mesh;

		//NonZeroCoefficients I_c_coeffs(coarsePb->HHO.nTotalCellUnknowns * coarsePb->HHO.nTotalFaceUnknowns); // Global interpolator from coarse faces onto coarse interiors
		NonZeroCoefficients J_f_c_coeffs(coarsePb->HHO.nTotalCellUnknowns * finePb->HHO.nTotalCellUnknowns); // Straight injector from coarse interior functions to fine interior functions
		//NonZeroCoefficients Pi_f_coeffs(finePb->HHO.nTotalCellUnknowns * finePb->HHO.nTotalFaceUnknowns); // Global projector from fine interiors onto fine faces

		int nCoarseFaceUnknowns = coarsePb->HHO.nLocalFaceUnknowns;
		int nCoarseCellUnknowns = coarsePb->HHO.nLocalCellUnknowns;

		int nFineFaceUnknowns = finePb->HHO.nLocalFaceUnknowns;
		int nFineCellUnknowns = finePb->HHO.nLocalCellUnknowns;

		for (auto ce : coarseMesh->Elements)
		{
			Poisson_HHO_Element<Dim>* coarseElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(ce);

			// I_c contribution
			/*Eigen::MatrixXd local_I_c = coarseElement->ComputeInterpolationMatrixFromFaces();
			for (auto face : coarseElement->Faces)
			{
				if (face->IsDomainBoundary)
					continue;

				BigNumber elemGlobalNumber = coarseElement->Number;
				BigNumber faceGlobalNumber = face->Number;
				BigNumber faceLocalNumber = coarseElement->LocalNumberOf(face);

				I_c_coeffs.Add(elemGlobalNumber*nCoarseCellUnknowns, faceGlobalNumber*nCoarseFaceUnknowns, local_I_c.block(0, faceLocalNumber*nCoarseFaceUnknowns, nCoarseCellUnknowns, nCoarseFaceUnknowns));
			}*/
			
			// J_f_c contribution
			Eigen::MatrixXd local_J_f_c = coarseElement->ComputeCanonicalInjectionMatrixCoarseToFine();
			for (auto fineElement : coarseElement->FinerElements)
			{
				BigNumber coarseElemGlobalNumber = coarseElement->Number;
				BigNumber fineElemGlobalNumber = fineElement->Number;
				BigNumber fineElemLocalNumber = coarseElement->LocalNumberOf(fineElement);

				J_f_c_coeffs.Add(fineElemGlobalNumber*nFineCellUnknowns, coarseElemGlobalNumber*nCoarseCellUnknowns, local_J_f_c.block(fineElemLocalNumber*nFineCellUnknowns, 0, nFineCellUnknowns, nCoarseCellUnknowns));
			}

			// Pi_f contribution
			/*for (auto fe : coarseElement->FinerElements)
			{
				Poisson_HHO_Element<Dim>* fineElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(fe);

				Eigen::MatrixXd local_Pi_f = fineElement->ComputeProjectorMatrixFromCellOntoFaces();
				for (auto face : fineElement->Faces)
				{
					if (face->IsDomainBoundary)
						continue;

					BigNumber elemGlobalNumber = fineElement->Number;
					BigNumber faceGlobalNumber = face->Number;
					BigNumber faceLocalNumber = fineElement->LocalNumberOf(face);

					Pi_f_coeffs.Add(faceGlobalNumber*nFineFaceUnknowns, elemGlobalNumber*nFineCellUnknowns, local_Pi_f.block(faceLocalNumber, 0, nFineFaceUnknowns, nFineCellUnknowns));
				}
			}*/
		}

		//Eigen::SparseMatrix<double> I_c(coarsePb->HHO.nTotalCellUnknowns, coarsePb->HHO.nTotalFaceUnknowns);
		Eigen::SparseMatrix<double> J_f_c(finePb->HHO.nTotalCellUnknowns, coarsePb->HHO.nTotalCellUnknowns);
		//Eigen::SparseMatrix<double> Pi_f(finePb->HHO.nTotalFaceUnknowns, finePb->HHO.nTotalCellUnknowns);

		//I_c_coeffs.Fill(I_c);
		J_f_c_coeffs.Fill(J_f_c);
		//Pi_f_coeffs.Fill(Pi_f);

		auto I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
		auto Pi_c = GetGlobalProjectorMatrixFromCellsToFaces(coarsePb);
		auto Pi_f = GetGlobalProjectorMatrixFromCellsToFaces(finePb);


		//cout << "------------------ " << " mult" << endl << (globalIc*globalPic) << endl;

		finePb->ExportMatrix(I_c, "I_c");
		finePb->ExportMatrix(J_f_c, "J_f_c");
		finePb->ExportMatrix(Pi_f, "Pi_f");

		finePb->ExportMatrix(Pi_c, "Pi_c");

		// Prolongation
		P = Pi_f * J_f_c * I_c;
		finePb->ExportMatrix(P, "P");
	}

	void SetupRestriction()
	{
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

		Eigen::SparseMatrix<double> M_c = GetGlobalMassMatrix(coarsePb); // Coarse mass matrix
		Eigen::SparseMatrix<double> M_f = GetGlobalMassMatrix(finePb); // Fine mass matrix

		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.analyzePattern(M_c);
		solver.factorize(M_c);
		Eigen::SparseMatrix<double> I(M_c.rows(), M_c.cols());
		I.setIdentity();
		auto invM_c = solver.solve(I);

		R = invM_c * P.transpose() * M_f;

		finePb->ExportMatrix(R, "R");
	}

	Eigen::SparseMatrix<double> GetGlobalMassMatrix(Poisson_HHO<Dim>* problem)
	{
		int faceLocalUnknowns = problem->HHO.nLocalFaceUnknowns;
		NonZeroCoefficients M_coeffs(problem->HHO.nElements * problem->HHO.nLocalFaceUnknowns * problem->HHO.nLocalFaceUnknowns);
		for (auto e : problem->_mesh->Elements)
		{
			Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);
			for (auto f : element->Faces)
			{
				if (f->IsDomainBoundary)
					continue;

				Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

				double weight = 1;// face->GetDiameter() / element->FrontierMeasure();
				Eigen::MatrixXd M_face = face->GetMassMatrix();

				BigNumber faceGlobalNumber = face->Number;
				BigNumber faceLocalNumber = element->LocalNumberOf(face);

				M_coeffs.Add(faceGlobalNumber*faceLocalUnknowns, faceGlobalNumber*faceLocalUnknowns, weight * M_face);
			}
		}
		Eigen::SparseMatrix<double> M(problem->HHO.nTotalFaceUnknowns, problem->HHO.nTotalFaceUnknowns);
		M_coeffs.Fill(M);
		return M;
	}

	Eigen::SparseMatrix<double> GetGlobalInterpolationMatrixFromFacesToCells(Poisson_HHO<Dim>* problem)
	{
		NonZeroCoefficients I_coeffs(problem->HHO.nTotalCellUnknowns * problem->HHO.nTotalFaceUnknowns);
		int nCellUnknowns = problem->HHO.nLocalCellUnknowns;
		int nFaceUnknowns = problem->HHO.nLocalFaceUnknowns;

		for (auto e : problem->_mesh->Elements)
		{
			Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);

			// I_c contribution
			Eigen::MatrixXd local_I = element->ComputeInterpolationMatrixFromFaces();

			//cout << "------------------ Element " << coarseElement->Number << " local_I_c" << endl << local_I_c << endl;
			for (auto face : element->Faces)
			{
				if (face->IsDomainBoundary)
					continue;

				BigNumber elemGlobalNumber = element->Number;
				BigNumber faceGlobalNumber = face->Number;
				BigNumber faceLocalNumber = element->LocalNumberOf(face);

				I_coeffs.Add(elemGlobalNumber*nCellUnknowns, faceGlobalNumber*nFaceUnknowns, local_I.block(0, faceLocalNumber*nFaceUnknowns, nCellUnknowns, nFaceUnknowns));
			}
		}
		Eigen::SparseMatrix<double> I(problem->HHO.nTotalCellUnknowns, problem->HHO.nTotalFaceUnknowns);
		I_coeffs.Fill(I);
		return I;
	}

	Eigen::SparseMatrix<double> GetGlobalProjectorMatrixFromCellsToFaces(Poisson_HHO<Dim>* problem)
	{
		NonZeroCoefficients Pi_coeffs(problem->HHO.nTotalCellUnknowns * problem->HHO.nTotalFaceUnknowns);
		int nCellUnknowns = problem->HHO.nLocalCellUnknowns;
		int nFaceUnknowns = problem->HHO.nLocalFaceUnknowns;

		for (auto e : problem->_mesh->Elements)
		{
			Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);

			Eigen::MatrixXd local_Pi = element->ComputeProjectorMatrixFromCellOntoFaces();
			for (auto face : element->Faces)
			{
				if (face->IsDomainBoundary)
					continue;

				BigNumber elemGlobalNumber = element->Number;
				BigNumber faceGlobalNumber = face->Number;
				BigNumber faceLocalNumber = element->LocalNumberOf(face);

				Pi_coeffs.Add(faceGlobalNumber*nFaceUnknowns, elemGlobalNumber*nCellUnknowns, local_Pi.block(faceLocalNumber*nFaceUnknowns, 0, nFaceUnknowns, nCellUnknowns));
			}
		}
		/*cout << "*****************************************************" << endl;
		for (auto face : problem->_mesh->Faces)
		{
			if (face->IsDomainBoundary)
				continue;

			Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(face->Element2);
			Eigen::MatrixXd local_Pi = element->ComputeProjectorMatrixFromCellOntoFaces();

			cout << "Face " << face->Number << " takes value from element " << element->Number << endl;

			BigNumber elemGlobalNumber = element->Number;
			BigNumber faceGlobalNumber = face->Number;
			BigNumber faceLocalNumber = element->LocalNumberOf(face);

			Pi_coeffs.Add(faceGlobalNumber*nFaceUnknowns, elemGlobalNumber*nCellUnknowns, local_Pi.block(faceLocalNumber*nFaceUnknowns, 0, nFaceUnknowns, nCellUnknowns));
		}*/

		Eigen::SparseMatrix<double> Pi(problem->HHO.nTotalFaceUnknowns, problem->HHO.nTotalCellUnknowns);
		Pi_coeffs.Fill(Pi);
		return Pi;
	}
};

template <int Dim>
class MultigridForHHO : public Multigrid
{
private:
	Poisson_HHO<Dim>* _problem;
	int _nLevels;
public:
	MultigridForHHO(Poisson_HHO<Dim>* problem, int nLevels) : Multigrid()
	{
		_problem = problem;
		_nLevels = nLevels;
	}

	void Setup(const Eigen::SparseMatrix<double>& A) override
	{
		IterativeSolver::Setup(A);

		LevelForHHO<Dim>* finerLevel = NULL;
		Poisson_HHO<Dim>* problem = _problem;
		for (int levelNumber = 0; levelNumber < _nLevels; levelNumber++)
		{
			LevelForHHO<Dim>* level = new LevelForHHO<Dim>(levelNumber, problem);
			if (levelNumber == 0)
				this->_fineLevel = level;
			else
			{
				finerLevel->CoarserLevel = level;
				level->FinerLevel = finerLevel;
			}

			finerLevel = level;

			if (levelNumber < _nLevels - 1)
			{
				dynamic_cast<CartesianGrid2D*>(problem->_mesh)->BuildCoarserMesh();
				problem = problem->GetProblemOnCoarserMesh();
				problem->Assemble(Action::None);
				/*for (auto f : problem->_mesh->Faces)
				{
					IntervalFace* face = dynamic_cast<IntervalFace*>(f);
					cout << *face << endl;
				}*/
			}
		}

		LevelForHHO<Dim>* level = dynamic_cast<LevelForHHO<Dim>*>(this->_fineLevel);
		while (level != NULL)
		{
			level->Setup();
			level = dynamic_cast<LevelForHHO<Dim>*>(level->CoarserLevel);
		}
	}

	/*GeometricLevel<Dim>* CreateGeometricLevel(int number, Mesh<Dim>* mesh)
	{
		return new LevelForHHO<Dim>(number, mesh);
	}*/
};