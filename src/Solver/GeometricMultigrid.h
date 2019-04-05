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
		Level* finerLevel = NULL;
		for (unsigned int i = 0; i < meshSequence.size(); i++)
		{
			Level* level = CreateLevel(i, meshSequence[i]);
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
	virtual GeometricLevel<Dim> CreateGeometricLevel(int number, Mesh<Dim> mesh) = 0;
};