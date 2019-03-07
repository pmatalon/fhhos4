#pragma once
#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "../Problem.h"
#include "../Mesh/Mesh.h"
#include "../Mesh/Face.h"
#include "Poisson_HHO_Element.h"
#include "Reconstructor.h"
#include "../Utils/NonZeroCoefficients.h"
#include "../Utils/L2.h"
using namespace std;

template <short Dim>
class Poisson_HHO : public Problem
{
private:
	Mesh<Dim>* _mesh;
	SourceFunction* _sourceFunction;
	FunctionalBasis<Dim>* _reconstructionBasis;
	FunctionalBasis<Dim>* _cellBasis;
	FunctionalBasis<Dim - 1>* _faceBasis;

public:
	Eigen::VectorXd ReconstructedSolution;

	Poisson_HHO(Mesh<Dim>* mesh, string solutionName, SourceFunction* sourceFunction, FunctionalBasis<Dim>* reconstructionBasis, FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim - 1>* faceBasis, string outputDirectory)
		: Problem(solutionName, outputDirectory)
	{	
		this->_mesh = mesh;
		this->_sourceFunction = sourceFunction;
		this->_reconstructionBasis = reconstructionBasis;
		this->_cellBasis = cellBasis;
		this->_faceBasis = faceBasis;
	}

	void Assemble(int penalizationCoefficient, Action action)
	{
		Mesh<Dim>* mesh = this->_mesh;
		FunctionalBasis<Dim>* reconstructionBasis = this->_reconstructionBasis;
		FunctionalBasis<Dim>* cellBasis = this->_cellBasis;
		FunctionalBasis<Dim - 1>* faceBasis = this->_faceBasis;


		cout << "Problem: Poisson " << Dim << "D" << endl;
		cout << "Subdivisions in each cartesian direction: " << mesh->N << endl;
		cout << "\tElements: " << mesh->Elements.size() << endl;
		cout << "\tFaces: " << mesh->Faces.size() << endl;
		cout << "Discretization: Hybrid High Order" << endl;
		cout << "\tReconstruction basis: " << reconstructionBasis->Name() << endl;
		cout << "\tCell basis: " << cellBasis->Name() << endl;
		cout << "\tFace basis: " << faceBasis->Name() << endl;

		BigNumber nTotalCellUnknowns = static_cast<int>(mesh->Elements.size()) * cellBasis->NumberOfLocalFunctionsInElement(NULL);
		BigNumber nTotalFaceUnknowns = static_cast<int>(mesh->Faces.size()) * faceBasis->NumberOfLocalFunctionsInElement(NULL);
		BigNumber nTotalHybridUnknowns = nTotalCellUnknowns + nTotalFaceUnknowns;
		cout << "Cell unknowns: " << nTotalCellUnknowns << " (" << cellBasis->NumberOfLocalFunctionsInElement(NULL) << " per cell)" << endl;
		cout << "Face unknowns: " << nTotalFaceUnknowns << " (" << faceBasis->NumberOfLocalFunctionsInElement(NULL) << " per face)" << endl;
		cout << "Total unknowns: " << nTotalHybridUnknowns << endl;

		bool autoPenalization = true;

		this->_fileName = "Poisson" + to_string(Dim) + "D" + this->_solutionName + "_n" + to_string(mesh->N) + "_HHO_" + reconstructionBasis->Name() + "_pen" + (autoPenalization ? "-1" : to_string(penalizationCoefficient));
		string matrixFilePath = this->_outputDirectory + "/" + this->_fileName + "_A.dat";
		string rhsFilePath = this->_outputDirectory + "/" + this->_fileName + "_b.dat";

		this->b = Eigen::VectorXd(nTotalHybridUnknowns);

		BigNumber nnzApproximate = mesh->Elements.size() * faceBasis->NumberOfLocalFunctionsInElement(NULL) * (2 * Dim + 1);
		NonZeroCoefficients consistencyCoeffs(nTotalHybridUnknowns, nTotalHybridUnknowns, nnzApproximate);
		NonZeroCoefficients stabilizationCoeffs(nTotalHybridUnknowns, nTotalHybridUnknowns, nnzApproximate);

		cout << "Assembly..." << endl;
		
		for (auto element : mesh->Elements)
		{
			//cout << "Element " << element->Number << endl;
			Poisson_HHO_Element<Dim>* hhoElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(element);
			hhoElement->InitReconstructor(reconstructionBasis, cellBasis, faceBasis);

			// Cell unknowns / Cell unknowns
			int nLocalCellUnknowns = cellBasis->NumberOfLocalFunctionsInElement(element);

			for (BasisFunction<Dim>* cellPhi1 : cellBasis->LocalFunctions)
			{
				BigNumber i = DOFNumber(element, cellPhi1);
				for (BasisFunction<Dim>* cellPhi2 : cellBasis->LocalFunctions)
				{
					BigNumber j = DOFNumber(element, cellPhi2);
					double consistencyTerm = hhoElement->ConsistencyTerm(cellPhi1, cellPhi2);
					consistencyCoeffs.Add(i, j, consistencyTerm);

					double stabilizationTerm = hhoElement->StabilizationTerm(cellPhi1, cellPhi2);
					stabilizationCoeffs.Add(i, j, stabilizationTerm);
				}
			}

			// Cell unknowns / Face unknowns
			for (BasisFunction<Dim>* cellPhi : cellBasis->LocalFunctions)
			{
				BigNumber i = DOFNumber(element, cellPhi);
				for (auto face : element->Faces)
				{
					for (BasisFunction<Dim-1>* facePhi : faceBasis->LocalFunctions)
					{
						BigNumber j = DOFNumber(face, facePhi);

						double consistencyTerm = hhoElement->ConsistencyTerm(face, cellPhi, facePhi);
						consistencyCoeffs.Add(i, j, consistencyTerm);
						consistencyCoeffs.Add(j, i, consistencyTerm);

						double stabilizationTerm = hhoElement->StabilizationTerm(face, cellPhi, facePhi);
						stabilizationCoeffs.Add(i, j, stabilizationTerm);
						stabilizationCoeffs.Add(j, i, stabilizationTerm);
					}
				}
			}

			// Face unknowns / Face unknowns
			for (auto face1 : element->Faces)
			{
				for (BasisFunction<Dim - 1>* facePhi1 : faceBasis->LocalFunctions)
				{
					BigNumber i = DOFNumber(face1, facePhi1);
					for (auto face2 : element->Faces)
					{
						for (BasisFunction<Dim - 1>* facePhi2 : faceBasis->LocalFunctions)
						{
							BigNumber j = DOFNumber(face2, facePhi2);

							double consistencyTerm = hhoElement->ConsistencyTerm(face1, facePhi1, face2, facePhi2);
							consistencyCoeffs.Add(i, j, consistencyTerm);
							//consistencyCoeffs.Add(j, i, consistencyTerm);

							double stabilizationTerm = hhoElement->StabilizationTerm(face1, facePhi1, face2, facePhi2);
							stabilizationCoeffs.Add(i, j, stabilizationTerm);
						}
					}
				}
			}

			// Right-hand side
			for (BasisFunction<Dim>* cellPhi : cellBasis->LocalFunctions)
			{
				BigNumber i = DOFNumber(element, cellPhi);
				this->b(i) = hhoElement->SourceTerm(cellPhi, this->_sourceFunction);
			}
			for (auto face : element->Faces)
			{
				for (BasisFunction<Dim - 1>* facePhi : faceBasis->LocalFunctions)
				{
					BigNumber i = DOFNumber(face, facePhi);
					this->b(i) = 0;
				}
			}

		}

		/*for (auto face : mesh->Faces)
		{
			if (face->IsDomainBoundary)
				continue;

			for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
			{
				BigNumber basisFunction1 = basis->GlobalFunctionNumber(face->Element1, phi1);
				for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
				{
				}
			}
		}*/

		Eigen::SparseMatrix<double> Acons = Eigen::SparseMatrix<double>(nTotalHybridUnknowns, nTotalHybridUnknowns);
		consistencyCoeffs.Fill(Acons);

		Eigen::SparseMatrix<double> Astab = Eigen::SparseMatrix<double>(nTotalHybridUnknowns, nTotalHybridUnknowns);
		stabilizationCoeffs.Fill(Astab);

		this->A = Acons + Astab;
		cout << "nnz(A) = " << this->A.nonZeros() << endl;

		cout << "Export..." << endl;
		Eigen::saveMarket(this->A, matrixFilePath);
		cout << "Matrix exported to \t" << matrixFilePath << endl;

		Eigen::saveMarketVector(this->b, rhsFilePath);
		cout << "RHS exported to \t" << rhsFilePath << endl;
	}

	void ReconstructSolution()
	{
		int nLocalReconstructUnknowns = this->_reconstructionBasis->NumberOfLocalFunctionsInElement(NULL);
		Eigen::VectorXd globalReconstructedSolution(this->_mesh->Elements.size() * nLocalReconstructUnknowns);

		BigNumber nTotalCellUnknowns = static_cast<int>(this->_mesh->Elements.size()) * this->_cellBasis->NumberOfLocalFunctionsInElement(NULL);

		for (auto element : this->_mesh->Elements)
		{
			Poisson_HHO_Element<Dim>* hhoElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(element);

			int nLocalCellUnknowns = this->_cellBasis->NumberOfLocalFunctionsInElement(element);
			int nLocalFaceUnknowns = this->_faceBasis->NumberOfLocalFunctionsInElement(NULL);

			Eigen::VectorXd localHybridSolution(nLocalCellUnknowns + nLocalFaceUnknowns * element->Faces.size());
			localHybridSolution.head(nLocalCellUnknowns) = this->Solution.segment(element->Number * nLocalCellUnknowns, nLocalCellUnknowns);
			for (auto face : element->Faces)
			{
				localHybridSolution.segment(nLocalCellUnknowns, nLocalFaceUnknowns) = this->Solution.segment(nTotalCellUnknowns + face->Number * nLocalFaceUnknowns, nLocalFaceUnknowns);
			}

			Eigen::VectorXd localReconstructedSolution = hhoElement->Reconstruct(localHybridSolution);

			globalReconstructedSolution.segment(element->Number * nLocalReconstructUnknowns, nLocalReconstructUnknowns) = localReconstructedSolution;
		}
		this->ReconstructedSolution = globalReconstructedSolution;
	}

	void ExtractSolution() override
	{
		Problem::ExtractSolution(this->ReconstructedSolution);
	}

private:
	BigNumber DOFNumber(Element<Dim>* element, BasisFunction<Dim>* cellPhi)
	{
		return element->Number * this->_cellBasis->Size() + cellPhi->LocalNumber;
	}
	BigNumber DOFNumber(Face<Dim>* face, BasisFunction<Dim-1>* facePhi)
	{
		BigNumber nTotalCellUnknowns = static_cast<int>(this->_mesh->Elements.size()) * this->_cellBasis->Size();
		return nTotalCellUnknowns + face->Number * this->_faceBasis->Size() + facePhi->LocalNumber;
	}
};

