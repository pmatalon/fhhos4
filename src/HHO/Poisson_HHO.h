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

public:

	Poisson_HHO(string solutionName) : Problem(solutionName)
	{	}

	void Assemble(Mesh<Dim>* mesh, FunctionalBasis<Dim>* elementBasis, FunctionalBasis<Dim-1>* faceBasis, string outputDirectory, Action action)
	{

		cout << "Discretization: Hybrid High Order" << endl;
		//cout << "\tBasis of polynomials: " << (dg->IsGlobalBasis() ? "global" : "") + basis->Name() << endl;

		/*cout << "Local functions: " << basis->NumberOfLocalFunctionsInElement(NULL) << endl;
		for (BasisFunction<Dim>* phi : basis->LocalFunctions)
			cout << "\t " << phi->ToString() << endl;*/
		BigNumber nElementUnknowns = static_cast<int>(mesh->Elements.size()) * elementBasis->NumberOfLocalFunctionsInElement(NULL);
		BigNumber nFaceUnknowns = static_cast<int>(mesh->Faces.size()) * faceBasis->NumberOfLocalFunctionsInElement(NULL);
		BigNumber nUnknowns = nElementUnknowns + nFaceUnknowns;
		cout << "Unknowns: " << nUnknowns << endl;

		string fileName = "Poisson" + to_string(Dim) + "D" + this->_solutionName + "_n" + to_string(mesh->N) + "_HHO_" + faceBasis->Name();
		string matrixFilePath = outputDirectory + "/" + fileName + "_A.dat";
		string rhsFilePath = outputDirectory + "/" + fileName + "_b.dat";

		this->b = Eigen::VectorXd(nUnknowns);

		BigNumber nnzApproximate = mesh->Elements.size() * faceBasis->NumberOfLocalFunctionsInElement(NULL) * (2 * Dim + 1);
		NonZeroCoefficients matrixCoeffs(nnzApproximate);

		cout << "Assembly..." << endl;

		//--------------------------------------------//
		// Iteration on the elements: diagonal blocks //
		//--------------------------------------------//

		for (auto element : mesh->Elements)
		{
			//cout << "Element " << element->Number << endl;
			Poisson_HHO_Element<Dim>* hhoElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(element);

			for (BasisFunction<Dim>* phi1 : faceBasis->LocalFunctions)
			{
				BigNumber basisFunction1 = faceBasis->GlobalFunctionNumber(element, phi1);

				// Current element (block diagonal)
				for (BasisFunction<Dim>* phi2 : faceBasis->LocalFunctions)
				{
					BigNumber basisFunction2 = faceBasis->GlobalFunctionNumber(element, phi2);

					//cout << "\t phi" << phi1->LocalNumber << " = " << phi1->ToString() << " phi" << phi2->LocalNumber << " = " << phi2->ToString() << endl;

					double interiorTerm = hhoElement->InteriorTerm(phi1, phi2);
					//cout << "\t\t interiorTerm = " << interiorTerm << endl;

					matrixCoeffs.Add(basisFunction1, basisFunction2, interiorTerm);
				}

				double rhs = hhoElement->SourceTerm(phi1);
				this->b(basisFunction1) = rhs;
			}
		}

		//---------------------------------------------//
		// Iteration on the faces: off-diagonal blocks //
		//---------------------------------------------//

		for (auto face : mesh->Faces)
		{
			if (face->IsDomainBoundary)
				continue;

			for (BasisFunction<Dim>* phi1 : faceBasis->LocalFunctions)
			{
				BigNumber basisFunction1 = faceBasis->GlobalFunctionNumber(face->Element1, phi1);
				for (BasisFunction<Dim>* phi2 : faceBasis->LocalFunctions)
				{
					BigNumber basisFunction2 = faceBasis->GlobalFunctionNumber(face->Element2, phi2);
					/*double coupling = dg->CouplingTerm(face, face->Element1, phi1, face->Element2, phi2);
					double penalization = dg->PenalizationTerm(face, face->Element1, phi1, face->Element2, phi2, penalizationCoefficient);

					if (extractMatrixComponents)
					{
						couplingCoeffs.Add(basisFunction1, basisFunction2, coupling);
						couplingCoeffs.Add(basisFunction2, basisFunction1, coupling);

						penCoeffs.Add(basisFunction1, basisFunction2, penalization);
						penCoeffs.Add(basisFunction2, basisFunction1, penalization);
					}
					matrixCoeffs.Add(basisFunction1, basisFunction2, coupling + penalization);
					matrixCoeffs.Add(basisFunction2, basisFunction1, coupling + penalization);*/
				}
			}
		}

		this->A = Eigen::SparseMatrix<double>(nUnknowns, nUnknowns);
		matrixCoeffs.Fill(this->A);
		cout << "nnz(A) = " << this->A.nonZeros() << endl;

		cout << "Export..." << endl;
		Eigen::saveMarket(this->A, matrixFilePath);
		cout << "Matrix exported to \t" << matrixFilePath << endl;

		Eigen::saveMarketVector(this->b, rhsFilePath);
		cout << "RHS exported to \t" << rhsFilePath << endl;

	}
};

