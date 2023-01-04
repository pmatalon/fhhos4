#pragma once
#include "../../FunctionalBasis/FunctionalBasis.h"

template <int Dim>
struct HHOParameters
{
	FunctionalBasis<Dim>* ReconstructionBasis;
	FunctionalBasis<Dim>* CellBasis;
	FunctionalBasis<Dim - 1>* FaceBasis;
	int OrthogonalizeElemBasesCode;
	int OrthogonalizeFaceBasesCode;

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
	BigNumber nTotalReconstructUnknowns;

	BigNumber nDirichletCoeffs;
	BigNumber nNeumannUnknowns;

	string Stabilization;

	HHOParameters(Mesh<Dim>* mesh, string stabilization, 
				FunctionalBasis<Dim>* reconstructionBasis, FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim - 1>* faceBasis, 
				int orthogonalizeElemBasesCode, int orthogonalizeFaceBasesCode)
	{
		this->ReconstructionBasis = reconstructionBasis;
		this->CellBasis = cellBasis;
		this->FaceBasis = faceBasis;
		this->OrthogonalizeElemBasesCode = orthogonalizeElemBasesCode;
		this->OrthogonalizeFaceBasesCode = orthogonalizeFaceBasesCode;

		nElements = mesh->Elements.size();
		nFaces = mesh->Faces.size();
		nInteriorFaces = mesh->InteriorFaces.size();
		nBoundaryFaces = mesh->BoundaryFaces.size();
		nDirichletFaces = mesh->DirichletFaces.size();
		nNeumannFaces = mesh->NeumannFaces.size();
		assert(nDirichletFaces + nNeumannFaces == nBoundaryFaces);
		assert(nBoundaryFaces + nInteriorFaces == nFaces);
		nInteriorAndNeumannFaces = nInteriorFaces + nNeumannFaces;

		nFaceUnknowns = faceBasis->Size();
		nCellUnknowns = cellBasis->Size();
		nReconstructUnknowns = reconstructionBasis->Size();

		nTotalCellUnknowns = nElements * nCellUnknowns;
		nTotalFaceUnknowns = nInteriorAndNeumannFaces * nFaceUnknowns;
		nTotalFaceCoeffs = nFaces * nFaceUnknowns;
		nTotalHybridUnknowns = nTotalCellUnknowns + nTotalFaceUnknowns;
		nTotalHybridCoeffs = nTotalCellUnknowns + nTotalFaceCoeffs;
		nTotalReconstructUnknowns = nElements * nReconstructUnknowns;

		nDirichletCoeffs = nDirichletFaces * nFaceUnknowns;
		nNeumannUnknowns = nNeumannFaces * nFaceUnknowns;

		Stabilization = stabilization;
	}

	bool OrthogonalizeElemBases()
	{
		return OrthogonalizeElemBasesCode > 0;
	}
	int NElemOrthogonalizations()
	{
		if (OrthogonalizeElemBasesCode == 1 || OrthogonalizeElemBasesCode == 3)
			return 1;
		else if (OrthogonalizeElemBasesCode == 2 || OrthogonalizeElemBasesCode == 4)
			return 2;
		return 0;
	}
	bool OrthonormalizeElemBases()
	{
		return OrthogonalizeElemBasesCode == 3 || OrthogonalizeElemBasesCode == 4;
	}

	bool OrthogonalizeFaceBases()
	{
		return OrthogonalizeFaceBasesCode > 0;
	}
	int NFaceOrthogonalizations()
	{
		if (OrthogonalizeFaceBasesCode == 1 || OrthogonalizeFaceBasesCode == 3)
			return 1;
		else if (OrthogonalizeFaceBasesCode == 2 || OrthogonalizeFaceBasesCode == 4)
			return 2;
		return 0;
	}
	bool OrthonormalizeFaceBases()
	{
		return OrthogonalizeFaceBasesCode == 3 || OrthogonalizeFaceBasesCode == 4;
	}
};