#pragma once
#include "../FunctionalBasis/FunctionalBasis.h"

template <int Dim>
struct HHOParameters
{
	FunctionalBasis<Dim>* ReconstructionBasis;
	FunctionalBasis<Dim>* CellBasis;
	FunctionalBasis<Dim - 1>* FaceBasis;
	int OrthogonalizeBasesCode;

	BigNumber nElements;
	BigNumber nFaces;
	BigNumber nInteriorFaces;
	BigNumber nBoundaryFaces;
	BigNumber nDirichletFaces;
	BigNumber nNeumannFaces;
	BigNumber nInteriorAndNeumannFaces;

	int nFaceUnknowns;
	int nCellUnknowns;
	int nReconstructUnknowns;

	BigNumber nTotalCellUnknowns;
	BigNumber nTotalFaceUnknowns;
	BigNumber nTotalFaceCoeffs;
	BigNumber nTotalHybridUnknowns;
	BigNumber nTotalHybridCoeffs;

	BigNumber nDirichletCoeffs;
	BigNumber nNeumannUnknowns;

	string Stabilization;

	HHOParameters(Mesh<Dim>* mesh, string stabilization, FunctionalBasis<Dim>* reconstructionBasis, FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim - 1>* faceBasis, int orthonormalizeBasesCode)
	{
		this->ReconstructionBasis = reconstructionBasis;
		this->CellBasis = cellBasis;
		this->FaceBasis = faceBasis;
		this->OrthogonalizeBasesCode = orthonormalizeBasesCode;

		nElements = mesh->Elements.size();
		nFaces = mesh->Faces.size();
		nInteriorFaces = mesh->InteriorFaces.size();
		nBoundaryFaces = mesh->BoundaryFaces.size();
		nDirichletFaces = mesh->DirichletFaces.size();
		nNeumannFaces = mesh->NeumannFaces.size();
		nInteriorAndNeumannFaces = nInteriorFaces + nNeumannFaces;

		nFaceUnknowns = faceBasis->Size();
		nCellUnknowns = cellBasis->Size();
		nReconstructUnknowns = reconstructionBasis->Size();

		nTotalCellUnknowns = nElements * nCellUnknowns;
		nTotalFaceUnknowns = nInteriorAndNeumannFaces * nFaceUnknowns;
		nTotalFaceCoeffs = nFaces * nFaceUnknowns;
		nTotalHybridUnknowns = nTotalCellUnknowns + nTotalFaceUnknowns;
		nTotalHybridCoeffs = nTotalCellUnknowns + nTotalFaceCoeffs;

		nDirichletCoeffs = nDirichletFaces * nFaceUnknowns;
		nNeumannUnknowns = nNeumannFaces * nFaceUnknowns;

		Stabilization = stabilization;
	}

	bool OrthogonalizeBases()
	{
		return OrthogonalizeBasesCode > 0;
	}
	int NOrthogonalizations()
	{
		if (OrthogonalizeBasesCode == 1 || OrthogonalizeBasesCode == 3)
			return 1;
		else if (OrthogonalizeBasesCode == 2 || OrthogonalizeBasesCode == 4)
			return 2;
		return 0;
	}
	bool OrthonormalizeBases()
	{
		return OrthogonalizeBasesCode == 3 || OrthogonalizeBasesCode == 4;
	}
};