#pragma once
#include "../PoissonProblem.h"
#include "Poisson_HHO_Element.h"
#include "../Utils/ElementParallelLoop.h"
using namespace std;

template <int Dim>
class Poisson_HHO : public PoissonProblem<Dim>
{
private:
	bool _staticCondensation = false;

	SparseMatrix _globalMatrix;
	Eigen::VectorXd _globalRHS;

public:
	HHOParameters<Dim>* HHO;
	Eigen::VectorXd ReconstructedSolution;

	Poisson_HHO(Mesh<Dim>* mesh, string rhsCode, SourceFunction* sourceFunction, HHOParameters<Dim>* hho, bool staticCondensation, DiffusionPartition<Dim>* diffusionPartition, string outputDirectory)
		: PoissonProblem<Dim>(mesh, diffusionPartition, rhsCode, sourceFunction, outputDirectory)
	{	
		this->HHO = hho;
		this->_staticCondensation = staticCondensation;

		this->_fileName = this->_fileName + "_HHO_" + HHO->ReconstructionBasis->Name() + (_staticCondensation ? "" : "_nostaticcond");

		// Re-numbering of the faces (interior first, then boundary)
		BigNumber faceNumber = 0;
		for (auto face : this->_mesh->InteriorFaces)
			face->Number = faceNumber++;
		for (auto face : this->_mesh->BoundaryFaces)
			face->Number = faceNumber++;
	}

	Poisson_HHO<Dim>* GetProblemOnCoarserMesh()
	{
		HHOParameters<Dim>* coarseHHO = new HHOParameters<Dim>(this->_mesh->CoarseMesh, HHO->Stabilization, HHO->ReconstructionBasis, HHO->CellBasis, HHO->FaceBasis);
		return new Poisson_HHO<Dim>(this->_mesh->CoarseMesh, this->_rhsCode, this->_sourceFunction, coarseHHO, _staticCondensation, this->_diffusionPartition, this->_outputDirectory);
	}

	void ExportFaces()
	{
		return ExportFaces("");
	}
	void ExportFaces(string suffix)
	{
		string filePath = this->GetFilePath("faces" + suffix);
		this->_mesh->ExportFacesToMatlab(filePath);
	}

	void PrintDiscretization()
	{
		cout << this->_mesh->Description() << endl;
		cout << "    Elements: " << HHO->nElements << endl;
		cout << "    Faces   : " << HHO->nFaces << " (" << HHO->nInteriorFaces << " interior + " << HHO->nBoundaryFaces << " boundary)" << endl;
		cout << "Discretization: Hybrid High Order (k = " << HHO->FaceBasis->GetDegree() << ")" << endl;
		cout << "    Reconstruction basis: " << HHO->ReconstructionBasis->Name() << endl;
		cout << "    Cell basis          : " << HHO->CellBasis->Name() << endl;
		cout << "    Face basis          : " << HHO->FaceBasis->Name() << endl;
		cout << "Cell unknowns : " << HHO->nTotalCellUnknowns << " (" << HHO->CellBasis->Size() << " per cell)" << endl;
		cout << "Face unknowns : " << HHO->nTotalFaceUnknowns << " (" << HHO->FaceBasis->Size() << " per interior face)" << endl;
		cout << "Total unknowns: " << HHO->nTotalHybridUnknowns << endl;
		cout << "System size   : " << (this->_staticCondensation ? HHO->nTotalFaceUnknowns : HHO->nTotalHybridUnknowns) << " (" << (this->_staticCondensation ? "statically condensed" : "no static condensation") << ")" << endl;
	}

	void Assemble(Action action)
	{
		Mesh<Dim>* mesh = this->_mesh;
		FunctionalBasis<Dim>* reconstructionBasis = HHO->ReconstructionBasis;
		FunctionalBasis<Dim>* cellBasis = HHO->CellBasis;
		FunctionalBasis<Dim - 1>* faceBasis = HHO->FaceBasis;

		if ((action & Action::LogAssembly) == Action::LogAssembly)
		{
			this->PrintPhysicalProblem();
			this->PrintDiscretization();
			cout << "Parallelism   : " << (BaseParallelLoop::GetDefaultNThreads() == 1 ? "sequential execution" : to_string(BaseParallelLoop::GetDefaultNThreads()) + " threads") << endl;
		}

		string matrixFilePath				= this->GetFilePath("A");
		string consistencyFilePath			= this->GetFilePath("A_cons");
		string stabilizationFilePath		= this->GetFilePath("A_stab");
		string reconstructionMatrixFilePath	= this->GetFilePath("Reconstruct");
		string rhsFilePath					= this->GetFilePath("b");

		if ((action & Action::LogAssembly) == Action::LogAssembly)
			cout << endl << "Assembly..." << endl;

		this->_globalRHS = Eigen::VectorXd(HHO->nTotalHybridUnknowns);
		this->_globalRHS.tail(HHO->nTotalFaceUnknowns) = Eigen::VectorXd::Zero(HHO->nTotalFaceUnknowns);

		// Compute some useful integrals on reference element and store them
		CartesianShape<Dim, Dim>::ReferenceShape.ComputeAndStoreCellMassMatrix(cellBasis);
		CartesianShape<Dim, Dim>::ReferenceShape.ComputeAndStoreReconstructMassMatrix(reconstructionBasis);
		CartesianShape<Dim, Dim>::ReferenceShape.ComputeAndStoreCellStiffnessMatrix(cellBasis);
		CartesianShape<Dim, Dim>::ReferenceShape.ComputeAndStoreReconstructK1StiffnessMatrix(this->_diffusionPartition->K1, reconstructionBasis);
		CartesianShape<Dim, Dim>::ReferenceShape.ComputeAndStoreReconstructK2StiffnessMatrix(this->_diffusionPartition->K2, reconstructionBasis);
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

			chunk->ThreadFuture = std::async([this, mesh, cellBasis, faceBasis, reconstructionBasis, action, chunk, &chunksMatrixCoeffs, &chunksConsistencyCoeffs, &chunksStabilizationCoeffs, &chunksReconstructionCoeffs]()
				{
					BigNumber nnzApproximate = chunk->Size() * (pow(HHO->nCellUnknowns, 2) + 4 * pow(HHO->nFaceUnknowns, 2) + HHO->nCellUnknowns * 4 * HHO->nFaceUnknowns);
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

						if (!this->_staticCondensation)
						{
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

		BigNumber nnzApproximate = mesh->Elements.size() * (pow(HHO->nCellUnknowns, 2) + 4 * pow(HHO->nFaceUnknowns, 2) + HHO->nCellUnknowns * 4 * HHO->nFaceUnknowns);
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

		SparseMatrix globalMatrix = SparseMatrix(HHO->nTotalHybridUnknowns, HHO->nTotalHybridUnknowns);
		matrixCoeffs.Fill(globalMatrix);

		this->_globalMatrix = globalMatrix;
		if (this->_staticCondensation)
		{
			if ((action & Action::LogAssembly) == Action::LogAssembly)
				cout << "Static condensation..." << endl;

			//SparseMatrix Att = this->_globalMatrix.topLeftCorner(HHO->nTotalCellUnknowns, HHO->nTotalCellUnknowns);
			SparseMatrix Aff = this->_globalMatrix.bottomRightCorner(HHO->nTotalFaceUnknowns, HHO->nTotalFaceUnknowns);
			SparseMatrix Atf = this->_globalMatrix.topRightCorner(HHO->nTotalCellUnknowns, HHO->nTotalFaceUnknowns);

			SparseMatrix inverseAtt = GetInverseAtt();

			this->A = Aff - Atf.transpose() * inverseAtt * Atf;

			Eigen::VectorXd bt = this->_globalRHS.head(HHO->nTotalCellUnknowns);
			Eigen::VectorXd bf = this->_globalRHS.tail(HHO->nTotalFaceUnknowns);
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
			SparseMatrix Acons = SparseMatrix(HHO->nTotalHybridUnknowns, HHO->nTotalHybridUnknowns);
			consistencyCoeffs.Fill(Acons);

			SparseMatrix Astab = SparseMatrix(HHO->nTotalHybridUnknowns, HHO->nTotalHybridUnknowns);
			stabilizationCoeffs.Fill(Astab);

			SparseMatrix reconstructionMatrix = SparseMatrix(HHO->nElements * reconstructionBasis->Size(), HHO->nTotalHybridCoeffs);
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
		// Init faces //
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->Faces, [this](Face<Dim>* f)
			{
				Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
				face->InitHHO(HHO);
			}
		);

		// Init Elements //
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this](Element<Dim>* e)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);
				element->InitHHO(HHO);
			}
		);
	}

	void ReconstructHigherOrderApproximation()
	{
		Eigen::VectorXd globalReconstructedSolution(HHO->nElements * HHO->nReconstructUnknowns);

		Eigen::VectorXd globalHybridSolution;

		if (this->_staticCondensation)
		{
			cout << "Solving cell unknowns..." << endl;
			Eigen::VectorXd facesSolution = this->Solution;

			SparseMatrix Atf = this->_globalMatrix.topRightCorner(HHO->nTotalCellUnknowns, HHO->nTotalFaceUnknowns);
			Eigen::VectorXd bt = this->_globalRHS.head(HHO->nTotalCellUnknowns);
			SparseMatrix inverseAtt = GetInverseAtt();

			globalHybridSolution = Eigen::VectorXd(HHO->nTotalHybridUnknowns);
			globalHybridSolution.tail(HHO->nTotalFaceUnknowns) = facesSolution;
			globalHybridSolution.head(HHO->nTotalCellUnknowns) = inverseAtt * (bt - Atf * facesSolution);
		}
		else
			globalHybridSolution = this->Solution;

		cout << "Reconstruction of higher order approximation..." << endl;
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &globalHybridSolution, &globalReconstructedSolution](Element<Dim>* e)
			{
				HHOParameters<Dim>* HHO = this->HHO;
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);

				Eigen::VectorXd localHybridSolution(HHO->nCellUnknowns + HHO->nFaceUnknowns * element->Faces.size());
				localHybridSolution.head(HHO->nCellUnknowns) = globalHybridSolution.segment(FirstDOFGlobalNumber(element), HHO->nCellUnknowns);
				for (auto face : element->Faces)
				{
					if (face->IsDomainBoundary)
						localHybridSolution.segment(element->FirstDOFNumber(face), HHO->nFaceUnknowns) = Eigen::VectorXd::Zero(HHO->nFaceUnknowns);
					else
						localHybridSolution.segment(element->FirstDOFNumber(face), HHO->nFaceUnknowns) = globalHybridSolution.segment(FirstDOFGlobalNumber(face), HHO->nFaceUnknowns);
				}

				Eigen::VectorXd localReconstructedSolution = element->Reconstruct(localHybridSolution);
				globalReconstructedSolution.segment(element->Number * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns) = localReconstructedSolution;
			});

		this->ReconstructedSolution = globalReconstructedSolution;
	}

	void ExtractTraceSystemSolution()
	{
		if (this->_staticCondensation)
		{
			Problem<Dim>::ExtractSolution(this->Solution, "Faces");
		}
	}

	void ExtractSolution() override
	{
		Problem<Dim>::ExtractSolution(this->ReconstructedSolution);
	}

	SparseMatrix GetInverseAtt()
	{
		ElementParallelLoop<Dim> parallelLoop(this->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nCellUnknowns * HHO->nCellUnknowns);
		parallelLoop.Execute([this](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);
				chunk->Results.Coeffs.Add(element->Number * HHO->nCellUnknowns, element->Number * HHO->nCellUnknowns, element->invAtt);
			});
		SparseMatrix inverseAtt = SparseMatrix(HHO->nTotalCellUnknowns, HHO->nTotalCellUnknowns);
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
		return element->Number * HHO->CellBasis->Size();
	}
	BigNumber FirstDOFGlobalNumber(Face<Dim>* face)
	{
		return this->HHO->nTotalCellUnknowns + face->Number * HHO->FaceBasis->Size();
	}
};

