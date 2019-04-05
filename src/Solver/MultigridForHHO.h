#pragma once
#include <Eigen/Sparse>
#include <vector>
#include "GeometricMultigrid.h"
#include "BlockSOR.h"
#include "../Mesh/Mesh.h"
using namespace std;
//using namespace Eigen;

template <int Dim>
class LevelForHHO : public GeometricLevel<Dim>
{
public:
	Eigen::SparseMatrix<double> R;
	Eigen::SparseMatrix<double> P;

	LevelForHHO(int number, Mesh<Dim>* mesh) 
		: GeometricLevel<Dim>(number, GaussSeidel(), ReverseGaussSeidel(), mesh)
	{}

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
		P = Eigen::SparseMatrix<double>(this->OperatorMatrixthis.rows(), this->CoarserLevel->OperatorMatrix.rows());

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

	GeometricLevel<Dim>* CreateGeometricLevel(int number, Mesh<Dim> mesh)
	{
		return new LevelForHHO<Dim>(number, mesh);
	}
};