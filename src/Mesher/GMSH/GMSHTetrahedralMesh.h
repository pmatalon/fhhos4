#pragma once
#include "GMSHMesh.h"
#include "../../Mesh/3D/TetrahedralMesh.h"
using namespace std;

class GMSHTetrahedralMesh : public GMSHMesh<3>, public TetrahedralMesh
{
public:
	GMSHTetrahedralMesh(TestCase<3>* testCase, string mshFile, BigNumber n = 0) :
		GMSHMesh(testCase, mshFile, n)
	{}
	GMSHTetrahedralMesh(TestCase<3>* testCase, string mshFile, string description, string fileNamePart, string geometryDescription, BigNumber n = 0) :
		GMSHMesh(testCase, mshFile, description, fileNamePart, geometryDescription, n)
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

	double AverageH() override
	{
		return GMSHMesh<3>::AverageH();
	}
protected:
	GMSHTetrahedralMesh(string description, string fileNamePart, string geometryDescription) 
		: GMSHMesh(description, fileNamePart)
	{
		this->_geometryDescription = geometryDescription;
	}

public:
	void RefineMeshBySplitting() override
	{
		GMSHMesh::RefineMeshBySplitting();
		this->FineMesh->ComesFrom.nFineElementsByCoarseElement = 8;
		this->FineMesh->ComesFrom.nFineFacesAddedByCoarseElement = 8;
		this->FineMesh->ComesFrom.nFineFacesByKeptCoarseFace = 4;
	}

	void RefineMesh(H_CoarsStgy strategy) override
	{
		if (this->FineMesh)
			assert(false && "Mesh already refined!");

		if (strategy == H_CoarsStgy::GMSHSplittingRefinement)
			this->RefineMeshBySplitting();
		else if (strategy == H_CoarsStgy::BeyRefinement)
			TetrahedralMesh::RefineMeshByBeyMethod();
		else
			GMSHMesh<3>::RefineMesh(strategy);
	}

protected:
	virtual GMSHMesh<3>* CreateNewGMSHMesh() override
	{
		return new GMSHTetrahedralMesh(this->_description, this->_fileNamePart, this->_geometryDescription);
	}
};