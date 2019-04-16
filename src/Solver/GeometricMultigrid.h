#pragma once
#include <Eigen/Sparse>
#include "Multigrid.h"
#include "Level.h"
#include "../Mesh/Mesh.h"
using namespace std;

template <int Dim>
class GeometricLevel : public Level
{
public:
	Mesh<Dim>* Grid;
	GeometricLevel(int number, Smoother preSmoother, Smoother postSmoother, Mesh<Dim>* mesh) : Level(number, preSmoother, postSmoother)
	{
		this->Grid = mesh;
	}
};

template <int Dim>
class GeometricMultigrid : public Multigrid
{
private:
	vector<Mesh<Dim>*> _meshSequence;
public:
	GeometricMultigrid(vector<Mesh<Dim>*> meshSequence) : Multigrid()
	{
		this->_meshSequence = meshSequence;
	}

	/*void Setup(const Eigen::SparseMatrix<double>& A) override
	{
		this->SetupLevelHierarchy();
		Multigrid::Setup(A);
	}*/
private:
	void SetupLevelHierarchy() override
	{
		Level* finerLevel = NULL;
		for (int i = 0; i < _meshSequence.size(); i++)
		{
			GeometricLevel<Dim>* level = this->CreateGeometricLevel(i, _meshSequence[i]);
			if (i == 0)
				this->_fineLevel = level;
			else
			{
				finerLevel->CoarserLevel = level;
				level->FinerLevel = finerLevel;
			}

			finerLevel = level;
		}
	}

protected:
	/*Level* CreateLevel(int number) override
	{
		return CreateGeometricLevel(number, _meshSequence[number]);
	}*/

	virtual GeometricLevel<Dim>* CreateGeometricLevel(int number, Mesh<Dim>* mesh) = 0;
};