#pragma once
#include "../GMSHMesh.h"
#include "TetrahedralMesh.h"
using namespace std;

class GMSHTetrahedralMesh : public GMSHMesh<3>, public TetrahedralMesh
{
public:
	GMSHTetrahedralMesh(string mshFile) :
		GMSHMesh(mshFile)
	{}
	GMSHTetrahedralMesh(string mshFile, string description, string fileNamePart) :
		GMSHMesh(mshFile, description, fileNamePart)
	{}

	string Description() override
	{
		return GMSHMesh<3>::Description();
	}

	string FileNamePart() override
	{
		return GMSHMesh<3>::FileNamePart();
	}

	double H() override
	{
		return GMSHMesh<3>::H();
	}

	double Regularity() override
	{
		return GMSHMesh<3>::Regularity();
	}
protected:
	GMSHTetrahedralMesh(string description, string fileNamePart) : GMSHMesh(description, fileNamePart)
	{}

public:
	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 8;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 8;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 4;
	}

	void RefineMesh(CoarseningStrategy strategy) override
	{
		if (this->FineMesh)
			assert(false && "Mesh already refined!");

		if (strategy == CoarseningStrategy::SplittingRefinement)
			this->RefineMeshBySplitting();
		else if (strategy == CoarseningStrategy::BeyRefinement)
			TetrahedralMesh::RefineMeshByBeyMethod();
		else
			GMSHMesh<3>::RefineMesh(strategy);
	}

protected:
	virtual GMSHMesh<3>* CreateEmptyGMSHMesh() override
	{
		return new GMSHTetrahedralMesh(this->_description, this->_fileNamePart);
	}
};