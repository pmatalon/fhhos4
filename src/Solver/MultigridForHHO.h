#pragma once
#include <Eigen/Sparse>
#include <vector>
#include "GeometricMultigrid.h"
#include "BlockSOR.h"
#include "../Mesh/Mesh.h"
#include "../HHO/Poisson_HHO_Element.h"
using namespace std;
//using namespace Eigen;

template <int Dim>
class LevelForHHO : public GeometricLevel<Dim>
{
private:
	Eigen::SparseMatrix<double> R;
	Eigen::SparseMatrix<double> P;

public:
	LevelForHHO(int number, Mesh<Dim>* mesh)
		: GeometricLevel<Dim>(number, Smoother(new GaussSeidel(), 1), Smoother(new ReverseGaussSeidel(), 1), mesh)
	{
		
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

private:
	void SetupRestriction() override
	{
		R = P.transpose();
	}

	void SetupProlongation() override
	{
		P = Eigen::SparseMatrix<double>(this->OperatorMatrix.rows(), this->CoarserLevel->OperatorMatrix.rows());
		Mesh<Dim>* coarseMesh = dynamic_cast<GeometricLevel<Dim>*>(this->CoarserLevel)->Grid;

		for (auto ce : coarseMesh->Elements)
		{
			Poisson_HHO_Element<Dim>* coarseElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(ce);
			Poisson_HHO_Element<Dim>* oneFinerElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(ce->FinerElements[0]);
			coarseElement->InitReconstructor(oneFinerElement->HHO()->ReconstructionBasis, oneFinerElement->HHO()->CellBasis, oneFinerElement->HHO()->FaceBasis); // TODO: make a light init (no nead to assemble the stabilization matrix for example)
			auto reconstructMatrix = coarseElement->HHO()->ReconstructionMatrix();
			auto reconstructionBasis = coarseElement->HHO()->ReconstructionBasis;
			auto faceBasis = coarseElement->HHO()->FaceBasis;

			for (auto f : coarseElement->FinerFacesRemoved)
			{
				/*// Passing to nodal values
				Eigen::MatrixXd modalToNodal_face(faceBasis->Size(), faceBasis->Size());
				vector<Point> facePoints = fineFace->GetNodalPoints(faceBasis);
				for (int i = 0; i < facePoints.size(); i++)
				{
					for (BasisFunction<Dim - 1>* phi : faceBasis->LocalFunctions)
						modalToNodal_face(i, phi->LocalNumber) = phi->Eval(facePoints[i]);
				}

				Eigen::MatrixXd modalToNodal_reconstruct(reconstructionBasis->Size(), reconstructionBasis->Size());
				vector<Point> reconstructPoints = coarseElement->GetNodalPoints(reconstructionBasis);
				for (int i = 0; i < reconstructPoints.size(); i++)
				{
					for (BasisFunction<Dim>* phi : reconstructionBasis->LocalFunctions)
						modalToNodal_reconstruct(i, phi->LocalNumber) = phi->Eval(reconstructPoints[i]);
				}

				reconstructedPolynomialRestrictedOnFace = modalToNodal_reconstruct * reconstructMatrix;*/
				Poisson_HHO_Face<Dim>* fineFace = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
				Eigen::MatrixXd ProjF = fineFace->GetProjFromReconstruct(this->_element);
			}
		}
	}

	void SetupOperator(const Eigen::SparseMatrix<double>& finerLevelOperatorMatrix) override
	{
		if (this->IsFinestLevel())
			this->OperatorMatrix = finerLevelOperatorMatrix;
		else
			this->OperatorMatrix = R * finerLevelOperatorMatrix * P;
	}
};

template <int Dim>
class MultigridForHHO : public GeometricMultigrid<Dim>
{
public:
	MultigridForHHO(vector<Mesh<Dim>*> meshSequence) : GeometricMultigrid<Dim>(meshSequence) {}

	GeometricLevel<Dim>* CreateGeometricLevel(int number, Mesh<Dim>* mesh)
	{
		return new LevelForHHO<Dim>(number, mesh);
	}
};