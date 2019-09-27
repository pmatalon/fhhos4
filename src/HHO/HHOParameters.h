#pragma once
#include "../FunctionalBasis/FunctionalBasis.h"

template <int Dim>
struct HHOParameters
{
	FunctionalBasis<Dim>* ReconstructionBasis;
	FunctionalBasis<Dim>* CellBasis;
	FunctionalBasis<Dim - 1>* FaceBasis;

	BigNumber nElements;
	BigNumber nFaces;
	BigNumber nInteriorFaces;
	BigNumber nBoundaryFaces;

	int nFaceUnknowns;
	int nCellUnknowns;
	int nReconstructUnknowns;

	BigNumber nTotalCellUnknowns;
	BigNumber nTotalFaceUnknowns;
	BigNumber nTotalHybridUnknowns;
	BigNumber nTotalHybridCoeffs;

	string Stabilization;

	HHOParameters(Mesh<Dim>* mesh, string stabilization, FunctionalBasis<Dim>* reconstructionBasis, FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim - 1>* faceBasis)
	{
		this->ReconstructionBasis = reconstructionBasis;
		this->CellBasis = cellBasis;
		this->FaceBasis = faceBasis;

		nElements = mesh->Elements.size();
		nFaces = mesh->Faces.size();
		nInteriorFaces = mesh->InteriorFaces.size();
		nBoundaryFaces = mesh->BoundaryFaces.size();

		nFaceUnknowns = faceBasis->Size();
		nCellUnknowns = cellBasis->Size();
		nReconstructUnknowns = reconstructionBasis->Size();

		nTotalCellUnknowns = nElements * nCellUnknowns;
		nTotalFaceUnknowns = nInteriorFaces * nFaceUnknowns;
		nTotalHybridUnknowns = nTotalCellUnknowns + nTotalFaceUnknowns;
		nTotalHybridCoeffs = nTotalCellUnknowns + nFaces * nFaceUnknowns;

		Stabilization = stabilization;
	}
};