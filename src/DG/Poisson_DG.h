#pragma once
#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "../Problem.h"
#include "../Mesh/Mesh.h"
#include "../Mesh/CartesianShape.h"
#include "../Utils/NonZeroCoefficients.h"
#include "../Utils/L2.h"
#include "../Utils/Action.h"
#include "../Utils/ParallelLoop.h"
#include "Poisson_DG_Element.h"
#include "Poisson_DG_Face.h"
using namespace std;

template <int Dim>
class Poisson_DG : public Problem
{
private:
	DiffusionPartition _diffusionPartition;
	SourceFunction* _sourceFunction;
public:

	Poisson_DG(string solutionName, SourceFunction* sourceFunction, DiffusionPartition diffusionPartition, string outputDirectory)
		: Problem(solutionName, outputDirectory), _diffusionPartition(diffusionPartition)
	{
		this->_sourceFunction = sourceFunction;
	}

	void Assemble(Mesh<Dim>* mesh, FunctionalBasis<Dim>* basis, int penalizationCoefficient, Action action)
	{
		cout << "Problem: Poisson " << Dim << "D";
		if (this->_diffusionPartition.Kappa1 != this->_diffusionPartition.Kappa2)
			cout << ", heterogeneous diffusion coefficient (k1=" << this->_diffusionPartition.Kappa1 << ", k2=" << this->_diffusionPartition.Kappa2 << ")";
		cout << endl;
		cout << "Analytical solution: " ;
		if (this->_solutionName.compare("sine") == 0)
			cout << "sine function";
		else if (this->_solutionName.compare("sine") == 0)
			cout << "polynomial function";
		else if (this->_solutionName.compare("hetero") == 0)
			cout << "heterogeneous-specific piecewise polynomial function";
		else
			cout << "unknown";
		cout << endl;
		cout << "Subdivisions in each cartesian direction: " << mesh->N << endl;
		cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
		cout << "\tPolynomial space: " << (basis->FullTensorization ? "Q" : "P") << endl;
		cout << "\tPolynomial basis: " << basis->Name() << endl;

		bool autoPenalization = penalizationCoefficient == -1;
		if (autoPenalization)
			penalizationCoefficient = pow(Dim, 2) * pow(basis->GetDegree() + 1, 2) * mesh->N; // Ralph-Hartmann
		cout << "\tPenalization coefficient: " << penalizationCoefficient << (autoPenalization ? " (automatic)" : "") << endl;
		
		cout << "Local functions: " << basis->NumberOfLocalFunctionsInElement(NULL) << endl;
		for (BasisFunction<Dim>* phi : basis->LocalFunctions)
			cout << "\t " << phi->ToString() << endl;
		BigNumber nUnknowns = static_cast<int>(mesh->Elements.size()) * basis->NumberOfLocalFunctionsInElement(NULL);
		cout << "Unknowns: " << nUnknowns << endl;

		string kappaString = "";
		if (this->_diffusionPartition.Kappa1 != 1)
		{
			char res[16];
			sprintf(res, "_kappa%g", this->_diffusionPartition.Kappa1);
			kappaString = res;
		}
		this->_fileName = "Poisson" + to_string(Dim) + "D" + this->_solutionName + kappaString + "_n" + to_string(mesh->N) + "_DG_SIPG_" + basis->Name() + "_pen" + (autoPenalization ? "-1" : to_string(penalizationCoefficient));
		string matrixFilePath			= this->_outputDirectory + "/" + this->_fileName + "_A.dat";
		string matrixVolumicFilePath	= this->_outputDirectory + "/" + this->_fileName + "_A_volumic.dat";
		string matrixCouplingFilePath	= this->_outputDirectory + "/" + this->_fileName + "_A_coupling.dat";
		string matrixPenFilePath		= this->_outputDirectory + "/" + this->_fileName + "_A_pen.dat";
		string massMatrixFilePath		= this->_outputDirectory + "/" + this->_fileName + "_Mass.dat";
		string rhsFilePath				= this->_outputDirectory + "/" + this->_fileName + "_b.dat";

		this->b = Eigen::VectorXd(nUnknowns);

		cout << "--------------------------------------------------------" << endl;
		cout << "Assembly..." << endl;

		CartesianShape<Dim, Dim>::ReferenceShape.ComputeMassMatrix(basis);
		CartesianShape<Dim, Dim>::ReferenceShape.ComputeStiffnessMatrix(basis);

		//--------------------------------------------//
		// Iteration on the elements: diagonal blocks //
		//--------------------------------------------//

		ParallelLoop parallelLoop(mesh->Elements.size());

		vector<NonZeroCoefficients> chunksMatrixCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksMassMatrixCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksVolumicCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksCouplingCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksPenCoeffs(parallelLoop.NThreads);
		
		for (unsigned int threadNumber = 0; threadNumber < parallelLoop.NThreads; threadNumber++)
		{
			ParallelChunk* chunk = parallelLoop.Chunks[threadNumber];

			chunk->ThreadFuture = std::async([this, mesh, basis, action, penalizationCoefficient, chunk, &chunksMatrixCoeffs, &chunksMassMatrixCoeffs, &chunksVolumicCoeffs, &chunksCouplingCoeffs, &chunksPenCoeffs]()
			{
				BigNumber nnzApproximate = chunk->Size() * basis->Size() * (2 * Dim + 1);
				NonZeroCoefficients matrixCoeffs(nnzApproximate);
				NonZeroCoefficients massMatrixCoeffs((action & Action::ExtractMassMatrix) == Action::ExtractMassMatrix ? nnzApproximate : 0);
				NonZeroCoefficients volumicCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);
				NonZeroCoefficients couplingCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);
				NonZeroCoefficients penCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);

				for (BigNumber iElem = chunk->Start; iElem < chunk->End; iElem++)
				{
					Poisson_DG_Element<Dim>* element = dynamic_cast<Poisson_DG_Element<Dim>*>(mesh->Elements[iElem]);
					//cout << "Element " << element->Number << endl;

					for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
					{
						BigNumber basisFunction1 = basis->GlobalFunctionNumber(element, phi1);

						// Current element (block diagonal)
						for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
						{
							BigNumber basisFunction2 = basis->GlobalFunctionNumber(element, phi2);

							//cout << "\t phi" << phi1->LocalNumber << " = " << phi1->ToString() << " phi" << phi2->LocalNumber << " = " << phi2->ToString() << endl;

							double volumicTerm = element->VolumicTerm(phi1, phi2, this->_diffusionPartition);
							//cout << "\t\t volumic = " << volumicTerm << endl;

							double coupling = 0;
							double penalization = 0;
							for (Face<Dim>* f : element->Faces)
							{
								Poisson_DG_Face<Dim>* face = dynamic_cast<Poisson_DG_Face<Dim>*>(f);

								double c = face->CouplingTerm(element, phi1, element, phi2, this->_diffusionPartition);
								double p = face->PenalizationTerm(element, phi1, element, phi2, penalizationCoefficient, this->_diffusionPartition);
								coupling += c;
								penalization += p;
								//cout << "\t\t " << face->ToString() << ":\t c=" << c << "\tp=" << p << endl;
							}

							//cout << "\t\t TOTAL = " << volumicTerm + coupling + penalization << endl;

							if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
							{
								volumicCoeffs.Add(basisFunction1, basisFunction2, volumicTerm);
								couplingCoeffs.Add(basisFunction1, basisFunction2, coupling);
								penCoeffs.Add(basisFunction1, basisFunction2, penalization);
							}
							matrixCoeffs.Add(basisFunction1, basisFunction2, volumicTerm + coupling + penalization);
							if ((action & Action::ExtractMassMatrix) == Action::ExtractMassMatrix)
							{
								double massTerm = element->MassTerm(phi1, phi2);
								massMatrixCoeffs.Add(basisFunction1, basisFunction2, massTerm);
							}
						}

						double rhs = element->SourceTerm(phi1, this->_sourceFunction);
						this->b(basisFunction1) = rhs;
					}
				}

				chunksMatrixCoeffs[chunk->ThreadNumber] = matrixCoeffs;
				chunksMassMatrixCoeffs[chunk->ThreadNumber] = massMatrixCoeffs;
				chunksVolumicCoeffs[chunk->ThreadNumber] = volumicCoeffs;
				chunksCouplingCoeffs[chunk->ThreadNumber] = couplingCoeffs;
				chunksPenCoeffs[chunk->ThreadNumber] = penCoeffs;
			}
			);
		}

		BigNumber nnzApproximate = mesh->Elements.size() * basis->Size() * (2 * Dim + 1);
		NonZeroCoefficients matrixCoeffs(nnzApproximate);
		NonZeroCoefficients massMatrixCoeffs((action & Action::ExtractMassMatrix) == Action::ExtractMassMatrix ? nnzApproximate : 0);
		NonZeroCoefficients volumicCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);
		NonZeroCoefficients couplingCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);
		NonZeroCoefficients penCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);

		parallelLoop.Wait();

		for (unsigned int threadNumber = 0; threadNumber < parallelLoop.NThreads; threadNumber++)
		{
			matrixCoeffs.Add(chunksMatrixCoeffs[threadNumber]);
			massMatrixCoeffs.Add(chunksMassMatrixCoeffs[threadNumber]);
			volumicCoeffs.Add(chunksVolumicCoeffs[threadNumber]);
			couplingCoeffs.Add(chunksCouplingCoeffs[threadNumber]);
			penCoeffs.Add(chunksPenCoeffs[threadNumber]);
		}

		chunksMatrixCoeffs.clear();
		chunksMassMatrixCoeffs.clear();
		chunksVolumicCoeffs.clear();
		chunksCouplingCoeffs.clear();
		chunksPenCoeffs.clear();

		//---------------------------------------------//
		// Iteration on the faces: off-diagonal blocks //
		//---------------------------------------------//

		ParallelLoop parallelLoopFaces(mesh->Faces.size());

		chunksMatrixCoeffs = vector<NonZeroCoefficients>(parallelLoopFaces.NThreads);
		chunksCouplingCoeffs = vector<NonZeroCoefficients>(parallelLoopFaces.NThreads);
		chunksPenCoeffs = vector<NonZeroCoefficients>(parallelLoopFaces.NThreads);


		for (unsigned int threadNumber = 0; threadNumber < parallelLoopFaces.NThreads; threadNumber++)
		{
			ParallelChunk* chunk = parallelLoopFaces.Chunks[threadNumber];

			chunk->ThreadFuture = std::async([this, mesh, basis, action, penalizationCoefficient, chunk, &chunksMatrixCoeffs, &chunksCouplingCoeffs, &chunksPenCoeffs]()
			{
				BigNumber nnzApproximate = chunk->Size() * basis->Size() * (2 * Dim + 1);
				NonZeroCoefficients matrixCoeffs(nnzApproximate);
				NonZeroCoefficients couplingCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);
				NonZeroCoefficients penCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);

				for (BigNumber iElem = chunk->Start; iElem < chunk->End; ++iElem)
				{
					Poisson_DG_Face<Dim>* face = dynamic_cast<Poisson_DG_Face<Dim>*>(mesh->Faces[iElem]);
					if (face->IsDomainBoundary)
						continue;

					//cout << "Face " << face->Number << endl;

					for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
					{
						BigNumber basisFunction1 = basis->GlobalFunctionNumber(face->Element1, phi1);
						for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
						{
							//cout << "\t phi" << phi1->LocalNumber << " = " << phi1->ToString() << " phi" << phi2->LocalNumber << " = " << phi2->ToString() << endl;

							BigNumber basisFunction2 = basis->GlobalFunctionNumber(face->Element2, phi2);
							double coupling = face->CouplingTerm(face->Element1, phi1, face->Element2, phi2, this->_diffusionPartition);
							double penalization = face->PenalizationTerm(face->Element1, phi1, face->Element2, phi2, penalizationCoefficient, this->_diffusionPartition);

							//cout << "\t\t\t c=" << coupling << "\tp=" << penalization << endl;

							if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
							{
								couplingCoeffs.Add(basisFunction1, basisFunction2, coupling);
								couplingCoeffs.Add(basisFunction2, basisFunction1, coupling);

								penCoeffs.Add(basisFunction1, basisFunction2, penalization);
								penCoeffs.Add(basisFunction2, basisFunction1, penalization);
							}
							matrixCoeffs.Add(basisFunction1, basisFunction2, coupling + penalization);
							matrixCoeffs.Add(basisFunction2, basisFunction1, coupling + penalization);
						}
					}
				}

				chunksMatrixCoeffs[chunk->ThreadNumber] = matrixCoeffs;
				chunksCouplingCoeffs[chunk->ThreadNumber] = couplingCoeffs;
				chunksPenCoeffs[chunk->ThreadNumber] = penCoeffs;
			}
			);
		}

		parallelLoopFaces.Wait();

		for (unsigned int threadNumber = 0; threadNumber < parallelLoopFaces.NThreads; threadNumber++)
		{
			matrixCoeffs.Add(chunksMatrixCoeffs[threadNumber]);
			couplingCoeffs.Add(chunksCouplingCoeffs[threadNumber]);
			penCoeffs.Add(chunksPenCoeffs[threadNumber]);
		}

		//---------------//
		// Matrix export //
		//---------------//

		this->A = Eigen::SparseMatrix<double>(nUnknowns, nUnknowns);
		matrixCoeffs.Fill(this->A);
		cout << "nnz(A) = " << this->A.nonZeros() << endl;

		if ((action & Action::ExtractSystem) == Action::ExtractSystem)
		{
			cout << "Export..." << endl;
			Eigen::saveMarket(this->A, matrixFilePath);
			cout << "Matrix exported to \t" << matrixFilePath << endl;

			Eigen::saveMarketVector(this->b, rhsFilePath);
			cout << "RHS exported to \t" << rhsFilePath << endl;
		}

		if ((action & Action::ExtractMassMatrix) == Action::ExtractMassMatrix)
		{
			Eigen::SparseMatrix<double> M(nUnknowns, nUnknowns);
			massMatrixCoeffs.Fill(M);
			Eigen::saveMarket(M, massMatrixFilePath);
			cout << "Mass matrix exported to \t" << massMatrixFilePath << endl;
		}

		if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
		{
			Eigen::SparseMatrix<double> V(nUnknowns, nUnknowns);
			volumicCoeffs.Fill(V);
			Eigen::saveMarket(V, matrixVolumicFilePath);
			cout << "Volumic part exported to \t" << matrixVolumicFilePath << endl;

			Eigen::SparseMatrix<double> C(nUnknowns, nUnknowns);
			couplingCoeffs.Fill(C);
			Eigen::saveMarket(C, matrixCouplingFilePath);
			cout << "Coupling part exported to \t" << matrixCouplingFilePath << endl;

			Eigen::SparseMatrix<double> P(nUnknowns, nUnknowns);
			penCoeffs.Fill(P);
			Eigen::saveMarket(P, matrixPenFilePath);
			cout << "Penalization part exported to \t" << matrixPenFilePath << endl;
		}

	}
};

