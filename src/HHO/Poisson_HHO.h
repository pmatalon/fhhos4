#pragma once
#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "../Problem.h"
#include "../Mesh/Mesh.h"
#include "../Mesh/Face.h"
#include "Poisson_HHO_Element.h"
#include "../Utils/NonZeroCoefficients.h"
#include "../Utils/L2.h"
#include "../Utils/ElementParallelLoop.h"
using namespace std;


template <int Dim>
struct HHOInfo
{
	BigNumber nElements;
	BigNumber nFaces;
	BigNumber nInteriorFaces;
	BigNumber nBoundaryFaces;

	int nLocalFaceUnknowns;
	int nLocalCellUnknowns;
	int nLocalReconstructUnknowns;

	BigNumber nTotalCellUnknowns;
	BigNumber nTotalFaceUnknowns;
	BigNumber nTotalHybridUnknowns;
	BigNumber nTotalHybridCoeffs;

	HHOInfo(Mesh<Dim>* mesh, FunctionalBasis<Dim>* reconstructionBasis, FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim - 1>* faceBasis)
	{
		nElements = mesh->Elements.size();
		nFaces = mesh->Faces.size();
		nInteriorFaces = mesh->InteriorFaces.size();
		nBoundaryFaces = mesh->BoundaryFaces.size();

		nLocalFaceUnknowns = faceBasis->Size();
		nLocalCellUnknowns = cellBasis->Size();
		nLocalReconstructUnknowns = reconstructionBasis->Size();

		nTotalCellUnknowns = nElements * nLocalCellUnknowns;
		nTotalFaceUnknowns = nInteriorFaces * nLocalFaceUnknowns;
		nTotalHybridUnknowns = nTotalCellUnknowns + nTotalFaceUnknowns;
		nTotalHybridCoeffs = nTotalCellUnknowns + nFaces * nLocalFaceUnknowns;;
	}
};

template <int Dim>
class Poisson_HHO : public Problem
{
private:
	SourceFunction* _sourceFunction;
	FunctionalBasis<Dim>* _reconstructionBasis;
	FunctionalBasis<Dim>* _cellBasis;
	FunctionalBasis<Dim - 1>* _faceBasis;
	bool _staticCondensation = false;

	Eigen::SparseMatrix<double> _globalMatrix;
	Eigen::VectorXd _globalRHS;

public:
	Mesh<Dim>* _mesh;
	HHOInfo<Dim> HHO;
	Eigen::VectorXd ReconstructedSolution;

	Poisson_HHO(Mesh<Dim>* mesh, string solutionName, SourceFunction* sourceFunction, FunctionalBasis<Dim>* reconstructionBasis, FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim - 1>* faceBasis, bool staticCondensation, string outputDirectory)
		: Problem(solutionName, outputDirectory), HHO(mesh, reconstructionBasis, cellBasis, faceBasis)
	{	
		this->_mesh = mesh;
		this->_sourceFunction = sourceFunction;
		this->_reconstructionBasis = reconstructionBasis;
		this->_cellBasis = cellBasis;
		this->_faceBasis = faceBasis;
		this->_staticCondensation = staticCondensation;
	}

	Poisson_HHO<Dim>* GetProblemOnCoarserMesh()
	{
		return new Poisson_HHO<Dim>(_mesh->CoarserMesh, _solutionName, _sourceFunction, _reconstructionBasis, _cellBasis, _faceBasis, _staticCondensation, _outputDirectory);
	}

	void Assemble(Action action)
	{
		Mesh<Dim>* mesh = this->_mesh;
		FunctionalBasis<Dim>* reconstructionBasis = this->_reconstructionBasis;
		FunctionalBasis<Dim>* cellBasis = this->_cellBasis;
		FunctionalBasis<Dim - 1>* faceBasis = this->_faceBasis;
		HHOInfo<Dim> hho = this->HHO;

		if ((action & Action::LogAssembly) == Action::LogAssembly)
		{
			cout << "Problem: Poisson " << Dim << "D" << endl;
			cout << "Subdivisions in each cartesian direction: " << mesh->N << endl;
			cout << "\tElements: " << hho.nElements << endl;
			cout << "\tFaces   : " << hho.nFaces << " (" << hho.nInteriorFaces << " interior + " << hho.nBoundaryFaces << " boundary)" << endl;
			cout << "Discretization: Hybrid High Order" << endl;
			cout << "\tReconstruction basis: " << reconstructionBasis->Name() << endl;
			cout << "\tCell basis          : " << cellBasis->Name() << endl;
			cout << "\tFace basis          : " << faceBasis->Name() << endl;
			cout << "Cell unknowns : " << hho.nTotalCellUnknowns << " (" << cellBasis->Size() << " per cell)" << endl;
			cout << "Face unknowns : " << hho.nTotalFaceUnknowns << " (" << faceBasis->Size() << " per interior face)" << endl;
			cout << "Total unknowns: " << hho.nTotalHybridUnknowns << endl;
			cout << "System size   : " << (this->_staticCondensation ? hho.nTotalFaceUnknowns : hho.nTotalHybridUnknowns) << " (" << (this->_staticCondensation ? "statically condensed" : "no static condensation") << ")" << endl;
			cout << "Parallelism   : " << (BaseParallelLoop::GetDefaultNThreads() == 1 ? "sequential execution" : to_string(BaseParallelLoop::GetDefaultNThreads()) + " threads") << endl;
		}
		this->_fileName = "Poisson" + to_string(Dim) + "D" + this->_solutionName + "_n" + to_string(mesh->N) + "_HHO_" + reconstructionBasis->Name() + "_pen-1" + (_staticCondensation ? "_staticcond" : "");
		string matrixFilePath				= this->_outputDirectory + "/" + this->_fileName + "_A.dat";
		string consistencyFilePath			= this->_outputDirectory + "/" + this->_fileName + "_A_cons.dat";
		string stabilizationFilePath		= this->_outputDirectory + "/" + this->_fileName + "_A_stab.dat";
		string reconstructionMatrixFilePath	= this->_outputDirectory + "/" + this->_fileName + "_Reconstruct.dat";
		string rhsFilePath					= this->_outputDirectory + "/" + this->_fileName + "_b.dat";

		if ((action & Action::LogAssembly) == Action::LogAssembly)
			cout << endl << "Assembly..." << endl;

		this->_globalRHS = Eigen::VectorXd(hho.nTotalHybridUnknowns);

		// Compute some useful integrals on reference element and store them
		CartesianShape<Dim, Dim>::ReferenceShape.ComputeAndStoreCellMassMatrix(cellBasis);
		CartesianShape<Dim, Dim>::ReferenceShape.ComputeAndStoreReconstructMassMatrix(reconstructionBasis);
		CartesianShape<Dim, Dim>::ReferenceShape.ComputeAndStoreCellStiffnessMatrix(cellBasis);
		CartesianShape<Dim, Dim>::ReferenceShape.ComputeAndStoreReconstructStiffnessMatrix(reconstructionBasis);
		CartesianShape<Dim, Dim>::ReferenceShape.ComputeAndStoreCellReconstructMassMatrix(cellBasis, reconstructionBasis);
		CartesianShape<Dim, Dim - 1>::ReferenceShape.ComputeAndStoreFaceMassMatrix(faceBasis);

		this->InitHHO();

		//-------------------------------//
		// Parallel loop on the elements //
		//-------------------------------//

		ParallelLoop<Element<Dim>*, EmptyResultChunk> parallelLoop(mesh->Elements);

		vector<NonZeroCoefficients> chunksMatrixCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksConsistencyCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? parallelLoop.NThreads : 0);
		vector<NonZeroCoefficients> chunksStabilizationCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? parallelLoop.NThreads : 0);
		vector<NonZeroCoefficients> chunksReconstructionCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? parallelLoop.NThreads : 0);
		
		for (unsigned int threadNumber = 0; threadNumber < parallelLoop.NThreads; threadNumber++)
		{
			ParallelChunk<EmptyResultChunk>* chunk = parallelLoop.Chunks[threadNumber];

			chunk->ThreadFuture = std::async([this, &hho, mesh, cellBasis, faceBasis, reconstructionBasis, action, chunk, &chunksMatrixCoeffs, &chunksConsistencyCoeffs, &chunksStabilizationCoeffs, &chunksReconstructionCoeffs]()
				{
					BigNumber nnzApproximate = chunk->Size() * (pow(hho.nLocalCellUnknowns, 2) + 4 * pow(hho.nLocalFaceUnknowns, 2) + hho.nLocalCellUnknowns * 4 * hho.nLocalFaceUnknowns);
					NonZeroCoefficients matrixCoeffs(nnzApproximate);
					NonZeroCoefficients consistencyCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);
					NonZeroCoefficients stabilizationCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);
					NonZeroCoefficients reconstructionCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);

					for (BigNumber iElem = chunk->Start; iElem < chunk->End; iElem++)
					{
						//--------------//
						//   Assembly   //
						//--------------//

						//cout << "Element " << element->Number << endl;
						Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(mesh->Elements[iElem]);

						// Cell unknowns / Cell unknowns
						for (BasisFunction<Dim>* cellPhi1 : cellBasis->LocalFunctions)
						{
							BigNumber i = DOFNumber(element, cellPhi1);
							for (BasisFunction<Dim>* cellPhi2 : cellBasis->LocalFunctions)
							{
								BigNumber j = DOFNumber(element, cellPhi2);

								double matrixTerm = element->MatrixTerm(cellPhi1, cellPhi2);
								matrixCoeffs.Add(i, j, matrixTerm);
								if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
								{
									double consistencyTerm = element->ConsistencyTerm(cellPhi1, cellPhi2);
									consistencyCoeffs.Add(i, j, consistencyTerm);

									double stabilizationTerm = element->StabilizationTerm(cellPhi1, cellPhi2);
									stabilizationCoeffs.Add(i, j, stabilizationTerm);
								}
							}
						}

						// Cell unknowns / Face unknowns
						for (BasisFunction<Dim>* cellPhi : cellBasis->LocalFunctions)
						{
							BigNumber i = DOFNumber(element, cellPhi);
							for (auto face : element->Faces)
							{
								if (face->IsDomainBoundary)
									continue;

								for (BasisFunction<Dim - 1>* facePhi : faceBasis->LocalFunctions)
								{
									BigNumber j = DOFNumber(face, facePhi);

									double matrixTerm = element->MatrixTerm(face, cellPhi, facePhi);
									matrixCoeffs.Add(i, j, matrixTerm);
									matrixCoeffs.Add(j, i, matrixTerm);
									if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
									{
										double consistencyTerm = element->ConsistencyTerm(face, cellPhi, facePhi);
										consistencyCoeffs.Add(i, j, consistencyTerm);
										consistencyCoeffs.Add(j, i, consistencyTerm);

										double stabilizationTerm = element->StabilizationTerm(face, cellPhi, facePhi);
										stabilizationCoeffs.Add(i, j, stabilizationTerm);
										stabilizationCoeffs.Add(j, i, stabilizationTerm);
									}
								}
							}
						}

						// Face unknowns / Face unknowns
						for (auto face1 : element->Faces)
						{
							if (face1->IsDomainBoundary)
								continue;

							for (BasisFunction<Dim - 1>* facePhi1 : faceBasis->LocalFunctions)
							{
								BigNumber i = DOFNumber(face1, facePhi1);
								for (auto face2 : element->Faces)
								{
									if (face2->IsDomainBoundary)
										continue;

									for (BasisFunction<Dim - 1>* facePhi2 : faceBasis->LocalFunctions)
									{
										BigNumber j = DOFNumber(face2, facePhi2);

										double matrixTerm = element->MatrixTerm(face1, facePhi1, face2, facePhi2);
										matrixCoeffs.Add(i, j, matrixTerm);
										if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
										{
											double consistencyTerm = element->ConsistencyTerm(face1, facePhi1, face2, facePhi2);
											consistencyCoeffs.Add(i, j, consistencyTerm);

											double stabilizationTerm = element->StabilizationTerm(face1, facePhi1, face2, facePhi2);
											stabilizationCoeffs.Add(i, j, stabilizationTerm);
										}
									}
								}
							}
						}

						// Right-hand side
						for (BasisFunction<Dim>* cellPhi : cellBasis->LocalFunctions)
						{
							BigNumber i = DOFNumber(element, cellPhi);
							this->_globalRHS(i) = element->SourceTerm(cellPhi, this->_sourceFunction);
						}
						for (auto face : element->Faces)
						{
							if (face->IsDomainBoundary)
								continue;

							for (BasisFunction<Dim - 1>* facePhi : faceBasis->LocalFunctions)
							{
								BigNumber i = DOFNumber(face, facePhi);
								this->_globalRHS(i) = 0;
							}
						}

						if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
						{
							// Global reconstruction matrix (for export)
							for (BasisFunction<Dim>* reconstructPhi : reconstructionBasis->LocalFunctions)
							{
								BigNumber i = element->Number * reconstructionBasis->Size() + reconstructPhi->LocalNumber;
								for (BasisFunction<Dim>* cellPhi : cellBasis->LocalFunctions)
								{
									BigNumber j = DOFNumber(element, cellPhi);
									reconstructionCoeffs.Add(i, j, element->ReconstructionTerm(reconstructPhi, cellPhi));
								}
								for (auto face : element->Faces)
								{
									for (BasisFunction<Dim - 1>* facePhi : faceBasis->LocalFunctions)
									{
										BigNumber j = DOFNumber(face, facePhi);
										reconstructionCoeffs.Add(i, j, element->ReconstructionTerm(reconstructPhi, face, facePhi));
									}
								}
							}
						}
					}

					chunksMatrixCoeffs[chunk->ThreadNumber] = matrixCoeffs;
					if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
					{
						chunksConsistencyCoeffs[chunk->ThreadNumber] = consistencyCoeffs;
						chunksStabilizationCoeffs[chunk->ThreadNumber] = stabilizationCoeffs;
						chunksReconstructionCoeffs[chunk->ThreadNumber] = reconstructionCoeffs;
					}
				}
			);
		}

		//------------------------------------//
		// Aggregation of the parallel chunks //
		//------------------------------------//

		BigNumber nnzApproximate = mesh->Elements.size() * (pow(hho.nLocalCellUnknowns, 2) + 4 * pow(hho.nLocalFaceUnknowns, 2) + hho.nLocalCellUnknowns * 4 * hho.nLocalFaceUnknowns);
		NonZeroCoefficients matrixCoeffs(nnzApproximate);
		NonZeroCoefficients consistencyCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);
		NonZeroCoefficients stabilizationCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);
		NonZeroCoefficients reconstructionCoeffs((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices ? nnzApproximate : 0);

		parallelLoop.Wait();

		for (unsigned int threadNumber = 0; threadNumber < parallelLoop.NThreads; threadNumber++)
		{
			matrixCoeffs.Add(chunksMatrixCoeffs[threadNumber]);
			if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
			{
				consistencyCoeffs.Add(chunksConsistencyCoeffs[threadNumber]);
				stabilizationCoeffs.Add(chunksStabilizationCoeffs[threadNumber]);
				reconstructionCoeffs.Add(chunksReconstructionCoeffs[threadNumber]);
			}
		}

		chunksMatrixCoeffs.clear();
		if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
		{
			chunksConsistencyCoeffs.clear();
			chunksStabilizationCoeffs.clear();
			chunksReconstructionCoeffs.clear();
		}

		//--------------------------//
		// Creation of the matrices //
		//--------------------------//

		Eigen::SparseMatrix<double> globalMatrix = Eigen::SparseMatrix<double>(hho.nTotalHybridUnknowns, hho.nTotalHybridUnknowns);
		matrixCoeffs.Fill(globalMatrix);

		this->_globalMatrix = globalMatrix;
		if (this->_staticCondensation)
		{
			if ((action & Action::LogAssembly) == Action::LogAssembly)
				cout << "Static condensation..." << endl;

			//Eigen::SparseMatrix<double> Att = this->_globalMatrix.topLeftCorner(hho.nTotalCellUnknowns, hho.nTotalCellUnknowns);
			Eigen::SparseMatrix<double> Aff = this->_globalMatrix.bottomRightCorner(hho.nTotalFaceUnknowns, hho.nTotalFaceUnknowns);
			Eigen::SparseMatrix<double> Atf = this->_globalMatrix.topRightCorner(hho.nTotalCellUnknowns, hho.nTotalFaceUnknowns);

			Eigen::SparseMatrix<double> inverseAtt = GetInverseAtt();

			this->A = Aff - Atf.transpose() * inverseAtt * Atf;

			Eigen::VectorXd bt = this->_globalRHS.head(hho.nTotalCellUnknowns);
			Eigen::VectorXd bf = this->_globalRHS.tail(hho.nTotalFaceUnknowns);
			this->b = bf - Atf.transpose() * inverseAtt * bt;
		}
		else
		{
			this->A = this->_globalMatrix;
			this->b = this->_globalRHS;
		}

		if ((action & Action::LogAssembly) == Action::LogAssembly)
			cout << Utils::MatrixInfo(this->A, "A") << endl;

		//------------------//
		//      Export      //
		//------------------//

		if ((action & Action::ExtractSystem) == Action::ExtractSystem)
		{
			cout << "Export:" << endl;
			Eigen::saveMarket(this->A, matrixFilePath);
			cout << "Matrix exported to                " << matrixFilePath << endl;

			Eigen::saveMarketVector(this->b, rhsFilePath);
			cout << "RHS exported to                   " << rhsFilePath << endl;
		}

		if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
		{
			Eigen::SparseMatrix<double> Acons = Eigen::SparseMatrix<double>(hho.nTotalHybridUnknowns, hho.nTotalHybridUnknowns);
			consistencyCoeffs.Fill(Acons);

			Eigen::SparseMatrix<double> Astab = Eigen::SparseMatrix<double>(hho.nTotalHybridUnknowns, hho.nTotalHybridUnknowns);
			stabilizationCoeffs.Fill(Astab);

			Eigen::SparseMatrix<double> reconstructionMatrix = Eigen::SparseMatrix<double>(hho.nElements * reconstructionBasis->Size(), hho.nTotalHybridCoeffs);
			reconstructionCoeffs.Fill(reconstructionMatrix);

			Eigen::saveMarket(Acons, consistencyFilePath);
			cout << "Consistency part exported to      " << consistencyFilePath << endl;

			Eigen::saveMarket(Astab, stabilizationFilePath);
			cout << "Stabilization part exported to    " << stabilizationFilePath << endl;

			Eigen::saveMarket(reconstructionMatrix, reconstructionMatrixFilePath);
			cout << "Reconstruction matrix exported to " << reconstructionMatrixFilePath << endl;
		}
	}

	void InitHHO()
	{
		// Global numbering of the faces (interior first, then boundary)
		BigNumber faceNumber = 0;
		for (auto face : this->_mesh->InteriorFaces)
			face->Number = faceNumber++;
		for (auto face : this->_mesh->BoundaryFaces)
			face->Number = faceNumber++;

		// Init faces //
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->Faces, [this](Face<Dim>* f)
			{
				Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
				face->InitHHO(this->_reconstructionBasis, this->_cellBasis, this->_faceBasis);
			}
		);

		// Init Elements //
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this](Element<Dim>* e)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);
				element->InitHHO(this->_reconstructionBasis, this->_cellBasis, this->_faceBasis);
			}
		);
	}

	void ReconstructSolution()
	{
		HHOInfo<Dim> hho = this->HHO;

		Eigen::VectorXd globalReconstructedSolution(hho.nElements * hho.nLocalReconstructUnknowns);

		Eigen::VectorXd globalHybridSolution;

		if (this->_staticCondensation)
		{
			cout << "Solving cell unknowns..." << endl;
			Eigen::VectorXd facesSolution = this->Solution;

			Eigen::SparseMatrix<double> Atf = this->_globalMatrix.topRightCorner(hho.nTotalCellUnknowns, hho.nTotalFaceUnknowns);
			Eigen::VectorXd bt = this->_globalRHS.head(hho.nTotalCellUnknowns);
			Eigen::SparseMatrix<double> inverseAtt = GetInverseAtt();

			globalHybridSolution = Eigen::VectorXd(hho.nTotalHybridUnknowns);
			globalHybridSolution.tail(hho.nTotalFaceUnknowns) = facesSolution;
			globalHybridSolution.head(hho.nTotalCellUnknowns) = inverseAtt * (bt - Atf * facesSolution);
		}
		else
			globalHybridSolution = this->Solution;

		cout << "Reconstruction of higher order approximation..." << endl;
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, hho, &globalHybridSolution, &globalReconstructedSolution](Element<Dim>* e)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);

				Eigen::VectorXd localHybridSolution(hho.nLocalCellUnknowns + hho.nLocalFaceUnknowns * element->Faces.size());
				localHybridSolution.head(hho.nLocalCellUnknowns) = globalHybridSolution.segment(FirstDOFGlobalNumber(element), hho.nLocalCellUnknowns);
				for (auto face : element->Faces)
				{
					if (face->IsDomainBoundary)
						localHybridSolution.segment(element->FirstDOFNumber(face), hho.nLocalFaceUnknowns) = Eigen::VectorXd::Zero(hho.nLocalFaceUnknowns);
					else
						localHybridSolution.segment(element->FirstDOFNumber(face), hho.nLocalFaceUnknowns) = globalHybridSolution.segment(FirstDOFGlobalNumber(face), hho.nLocalFaceUnknowns);
				}

				Eigen::VectorXd localReconstructedSolution = element->Reconstruct(localHybridSolution);
				globalReconstructedSolution.segment(element->Number * hho.nLocalReconstructUnknowns, hho.nLocalReconstructUnknowns) = localReconstructedSolution;
			});

		this->ReconstructedSolution = globalReconstructedSolution;
	}

	void ExtractTraceSystemSolution()
	{
		if (this->_staticCondensation)
		{
			Problem::ExtractSolution(this->Solution, "Faces");
		}
	}

	void ExtractSolution() override
	{
		Problem::ExtractSolution(this->ReconstructedSolution);
	}

	/*void Solve() override
	{
		vector<Mesh<Dim>*> meshSequence(2);
		dynamic_cast<CartesianGrid2D*>(this->_mesh)->BuildCoarserMesh();
		meshSequence.push_back(this->_mesh);
		meshSequence.push_back(this->_mesh->CoarserMesh);
		MultigridForHHO<Dim> mg(meshSequence);
	}*/

	Eigen::SparseMatrix<double> GetInverseAtt()
	{
		ElementParallelLoop<Dim> parallelLoop(_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(this->HHO.nLocalCellUnknowns * this->HHO.nLocalCellUnknowns);
		parallelLoop.Execute([this](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);
				int nCellUnknowns = this->HHO.nLocalCellUnknowns;
				chunk->Results.Coeffs.Add(element->Number * nCellUnknowns, element->Number * nCellUnknowns, element->invAtt);
			});
		Eigen::SparseMatrix<double> inverseAtt = Eigen::SparseMatrix<double>(this->HHO.nTotalCellUnknowns, this->HHO.nTotalCellUnknowns);
		parallelLoop.Fill(inverseAtt);
		return inverseAtt;
	}

private:
	BigNumber DOFNumber(Element<Dim>* element, BasisFunction<Dim>* cellPhi)
	{
		return FirstDOFGlobalNumber(element) + cellPhi->LocalNumber;
	}
	BigNumber DOFNumber(Face<Dim>* face, BasisFunction<Dim-1>* facePhi)
	{
		return FirstDOFGlobalNumber(face) + facePhi->LocalNumber;
	}
	BigNumber FirstDOFGlobalNumber(Element<Dim>* element)
	{
		return element->Number * this->_cellBasis->Size();
	}
	BigNumber FirstDOFGlobalNumber(Face<Dim>* face)
	{
		return this->HHO.nTotalCellUnknowns + face->Number * this->_faceBasis->Size();
	}
};

