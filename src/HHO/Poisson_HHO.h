#pragma once
#include "../Problem/PoissonProblem.h"
#include "Poisson_HHO_Element.h"
#include "../Utils/ElementParallelLoop.h"
#include "../Mesh/2D/Triangle.h"
#include "../Mesh/3D/Tetrahedron.h"
using namespace std;

template <int Dim>
class Poisson_HHO : public PoissonProblem<Dim>
{
private:
	bool _staticCondensation = false;

	SparseMatrix _globalSystemMatrix;
	Vector _globalRHS;
	Vector _dirichletCond;
public:
	HHOParameters<Dim>* HHO;
	Vector ReconstructedSolution;
	Vector GlobalHybridSolution;

	Poisson_HHO(Mesh<Dim>* mesh, string rhsCode, SourceFunction* sourceFunction, HHOParameters<Dim>* hho, bool staticCondensation, DiffusionPartition<Dim>* diffusionPartition, BoundaryConditions* bc, string outputDirectory)
		: PoissonProblem<Dim>(mesh, diffusionPartition, rhsCode, sourceFunction, bc, outputDirectory)
	{	
		this->HHO = hho;
		this->_staticCondensation = staticCondensation;

		this->_fileName = this->_fileName + "_HHO_" + HHO->ReconstructionBasis->Name() + (_staticCondensation ? "" : "_nostaticcond");

		// Re-numbering of the faces: interior first, then Neumann, and Dirichlet at the end (because they will be eliminated from the system)
		BigNumber faceNumber = 0;
		for (auto face : this->_mesh->InteriorFaces)
			face->Number = faceNumber++;
		for (auto face : this->_mesh->NeumannFaces)
			face->Number = faceNumber++;
		for (auto face : this->_mesh->DirichletFaces)
			face->Number = faceNumber++;
	}

	Poisson_HHO<Dim>* GetProblemOnCoarserMesh()
	{
		HHOParameters<Dim>* coarseHHO = new HHOParameters<Dim>(this->_mesh->CoarseMesh, HHO->Stabilization, HHO->ReconstructionBasis, HHO->CellBasis, HHO->FaceBasis);
		return new Poisson_HHO<Dim>(this->_mesh->CoarseMesh, this->_rhsCode, this->_sourceFunction, coarseHHO, _staticCondensation, this->_diffusionPartition, this->_boundaryConditions, this->_outputDirectory);
	}

	double L2Error(DomFunction exactSolution) override
	{
		return Problem<Dim>::L2Error(HHO->ReconstructionBasis, this->ReconstructedSolution, exactSolution);
	}

	void AssertSchemeConvergence(double l2Error)
	{
		if (Dim == 1)
			return;

		double h = this->_mesh->H();
		int k = this->HHO->FaceBasis->GetDegree();

		double constant = ConvergenceHiddenConstant(k);
		double upperBound;
		if ((this->_rhsCode.compare("sine") == 0 || this->_rhsCode.compare("poly") == 0) && this->_diffusionPartition->IsHomogeneous) // solution is H^{r+2} with r in <= k
		{
			double r = k;
			if (k == 0)
				upperBound = constant * pow(h, 2);
			else
				upperBound = constant * pow(h, r + 2);
		}
		else if (this->_rhsCode.compare("kellogg") == 0) // solution is H^{1+eps}
		{
			double eps = 0.1;
			upperBound = constant * pow(h, 2 * eps); // source: Di Pietro's mail, but he's not totally sure...
		}
		else
			return;

		cout << "UpperBound = " << upperBound << endl;
		if (l2Error > upperBound)
			cout << Utils::BeginRed << "The L2 error is bigger than the theoretical upper bound." << Utils::EndColor << endl;
	}
private:
	double ConvergenceHiddenConstant(int k)
	{
		// Depends on the domain. Very rough numerical estimate for [0, 1]^Dim:
		//return pow(10, Dim);
		return 100;
	}

public:
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
		cout << "Mesh: " << this->_mesh->Description() << endl;
		cout << "    Elements: " << HHO->nElements << endl;
		cout << "    Faces   : " << HHO->nFaces << " (" << HHO->nInteriorFaces << " interior + " << HHO->nBoundaryFaces << " boundary)" << endl;
		cout << "    h       : " << scientific << this->_mesh->H() << defaultfloat << endl;
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

		this->_globalRHS = Vector(HHO->nTotalHybridUnknowns);
		this->_globalRHS.tail(HHO->nTotalFaceUnknowns) = Vector::Zero(HHO->nTotalFaceUnknowns);

		// Compute some useful integrals on reference element and store them
		// - Cartesian element
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreCellMassMatrix(cellBasis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreReconstructMassMatrix(reconstructionBasis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreCellStiffnessMatrix(cellBasis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreReconstructK1StiffnessMatrix(this->_diffusionPartition->K1, reconstructionBasis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreReconstructK2StiffnessMatrix(this->_diffusionPartition->K2, reconstructionBasis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix(cellBasis, reconstructionBasis);
		CartesianShape<Dim, Dim - 1>::InitReferenceShape()->ComputeAndStoreFaceMassMatrix(faceBasis);
		if (Dim == 2)
		{
			// - Triangle
			TriangleShape::InitReferenceShape()->ComputeAndStoreCellMassMatrix((FunctionalBasis<2>*)cellBasis);
			TriangleShape::InitReferenceShape()->ComputeAndStoreReconstructMassMatrix((FunctionalBasis<2>*)reconstructionBasis);
			//Triangle::RefTriangle.ComputeAndStoreCellStiffnessMatrix((FunctionalBasis<2>*)cellBasis);
			//Triangle::RefTriangle.ComputeAndStoreReconstructK1StiffnessMatrix(this->_diffusionPartition->K1, reconstructionBasis);
			//Triangle::RefTriangle.ComputeAndStoreReconstructK2StiffnessMatrix(this->_diffusionPartition->K2, reconstructionBasis);
			TriangleShape::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix((FunctionalBasis<2>*)cellBasis, (FunctionalBasis<2>*)reconstructionBasis);
		}
		else if (Dim == 3)
		{
			// - Tetrahedron
			TetrahedronShape::InitReferenceShape()->ComputeAndStoreCellMassMatrix((FunctionalBasis<3>*)cellBasis);
			TetrahedronShape::InitReferenceShape()->ComputeAndStoreReconstructMassMatrix((FunctionalBasis<3>*)reconstructionBasis);
			TetrahedronShape::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix((FunctionalBasis<3>*)cellBasis, (FunctionalBasis<3>*)reconstructionBasis);
		}

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
							//-------------//
							// Cell / Cell //
							//-------------//

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

						//-------------//
						// Cell / Face //
						//-------------//

						for (BasisFunction<Dim>* cellPhi : cellBasis->LocalFunctions)
						{
							BigNumber i = DOFNumber(element, cellPhi);
							for (auto face : element->Faces)
							{
								for (BasisFunction<Dim - 1>* facePhi : faceBasis->LocalFunctions)
								{
									BigNumber j = DOFNumber(face, facePhi);

									double matrixTerm = element->MatrixTerm(face, cellPhi, facePhi);
									matrixCoeffs.Add(i, j, matrixTerm);
									matrixCoeffs.Add(j, i, matrixTerm);
									if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices && !face->HasDirichletBC())
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

						//-------------//
						// Face / Face //
						//-------------//

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

										double matrixTerm = element->MatrixTerm(face1, facePhi1, face2, facePhi2);
										matrixCoeffs.Add(i, j, matrixTerm);
										if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices && !face1->HasDirichletBC() && !face2->HasDirichletBC())
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

						//-----------------//
						// Right-hand side //
						//-----------------//

						for (BasisFunction<Dim>* cellPhi : cellBasis->LocalFunctions)
						{
							BigNumber i = DOFNumber(element, cellPhi);
							this->_globalRHS(i) = element->SourceTerm(cellPhi, this->_sourceFunction);
						}

						//-------------------------------------------//
						// Global reconstruction matrix (for export) //
						//-------------------------------------------//

						if ((action & Action::ExtractComponentMatrices) == Action::ExtractComponentMatrices)
						{
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


		SparseMatrix extendedMatrix = SparseMatrix(HHO->nTotalHybridCoeffs, HHO->nTotalHybridCoeffs);
		matrixCoeffs.Fill(extendedMatrix);

		// Extended matrix where the Dirichlet "unknowns" are eliminated
		this->_globalSystemMatrix = extendedMatrix.topLeftCorner(HHO->nTotalHybridUnknowns, HHO->nTotalHybridUnknowns);

		//---------------------------------//
		// Boundary conditions enforcement //
		//---------------------------------//

		// Neumann
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->NeumannFaces, [this, faceBasis](Face<Dim>* f)
			{
				Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
				BigNumber i = FirstDOFGlobalNumber(f);
				this->_globalRHS.segment(i, HHO->nFaceUnknowns) = face->ProjectOnBasis(faceBasis, this->_boundaryConditions->NeumannFunction);
			}
		);

		// Dirichlet
		this->_dirichletCond = Vector(HHO->nDirichletUnknowns);
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->DirichletFaces, [this, faceBasis](Face<Dim>* f)
			{
				Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
				BigNumber i = FirstDOFGlobalNumber(f) - HHO->nTotalHybridUnknowns;
				this->_dirichletCond.segment(i, HHO->nFaceUnknowns) = face->InvFaceMassMatrix()*face->ProjectOnBasis(faceBasis, this->_boundaryConditions->DirichletFunction);
			}
		);
		//cout << this->_dirichletCond << endl;
		this->_globalRHS -= extendedMatrix.topRightCorner(HHO->nTotalHybridUnknowns, HHO->nDirichletUnknowns) * this->_dirichletCond;

		//---------------------//
		// Static condensation //
		//---------------------//

		if (this->_staticCondensation)
		{
			if ((action & Action::LogAssembly) == Action::LogAssembly)
				cout << "Static condensation..." << endl;

			SparseMatrix Aff = this->_globalSystemMatrix.bottomRightCorner(HHO->nTotalFaceUnknowns, HHO->nTotalFaceUnknowns);
			SparseMatrix Atf = this->_globalSystemMatrix.topRightCorner(HHO->nTotalCellUnknowns, HHO->nTotalFaceUnknowns);

			SparseMatrix inverseAtt = GetInverseAtt();

			this->A = (Aff - Atf.transpose() * inverseAtt * Atf).pruned();

			Vector bt = this->_globalRHS.head(HHO->nTotalCellUnknowns);
			Vector bf = this->_globalRHS.tail(HHO->nTotalFaceUnknowns);
			this->b = bf - Atf.transpose() * inverseAtt * bt;
		}
		else
		{
			this->A = this->_globalSystemMatrix;
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
		Vector globalReconstructedSolution(HHO->nElements * HHO->nReconstructUnknowns);
		this->GlobalHybridSolution = Vector(HHO->nTotalHybridCoeffs);

		if (this->_staticCondensation)
		{
			cout << "Solving cell unknowns..." << endl;
			Vector facesSolution = this->SystemSolution;

			SparseMatrix Atf = this->_globalSystemMatrix.topRightCorner(HHO->nTotalCellUnknowns, HHO->nTotalFaceUnknowns);
			Vector bt = this->_globalRHS.head(HHO->nTotalCellUnknowns);
			SparseMatrix inverseAtt = GetInverseAtt();

			this->GlobalHybridSolution.head(HHO->nTotalCellUnknowns) = inverseAtt * (bt - Atf * facesSolution);
			this->GlobalHybridSolution.segment(HHO->nTotalCellUnknowns, HHO->nTotalFaceUnknowns) = facesSolution;
		}
		else
			this->GlobalHybridSolution.head(HHO->nTotalHybridUnknowns) = this->SystemSolution;

		// Dirichlet boundary conditions
		this->GlobalHybridSolution.tail(HHO->nDirichletUnknowns) = this->_dirichletCond;



		cout << "Reconstruction of higher order approximation..." << endl;
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &globalReconstructedSolution](Element<Dim>* e)
			{
				HHOParameters<Dim>* HHO = this->HHO;
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);

				Vector localHybridSolution(HHO->nCellUnknowns + HHO->nFaceUnknowns * element->Faces.size());
				localHybridSolution.head(HHO->nCellUnknowns) = this->GlobalHybridSolution.segment(FirstDOFGlobalNumber(element), HHO->nCellUnknowns);
				for (auto face : element->Faces)
				{
					localHybridSolution.segment(element->FirstDOFNumber(face), HHO->nFaceUnknowns) = this->GlobalHybridSolution.segment(FirstDOFGlobalNumber(face), HHO->nFaceUnknowns);
				}

				Vector localReconstructedSolution = element->Reconstruct(localHybridSolution);
				globalReconstructedSolution.segment(element->Number * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns) = localReconstructedSolution;
			});

		this->ReconstructedSolution = globalReconstructedSolution;
	}

	void ExtractTraceSystemSolution()
	{
		if (this->_staticCondensation)
		{
			Problem<Dim>::ExtractSolution(this->SystemSolution, "Faces");
		}
	}

	void ExtractHybridSolution()
	{
		Problem<Dim>::ExtractSolution(this->GlobalHybridSolution, "Hybrid");
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

	Vector ProjectOnFaceDiscreteSpace(DomFunction func)
	{
		Vector vectorOfDoFs = Vector(HHO->nTotalFaceUnknowns);
		auto faceBasis = HHO->FaceBasis;
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->Faces, [this, faceBasis, &vectorOfDoFs, func](Face<Dim>* f)
			{
				if (f->HasDirichletBC())
					return;
				Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
				BigNumber i = face->Number * faceBasis->Size();// FirstDOFGlobalNumber(f);

				DenseMatrix m = face->InvFaceMassMatrix();
				Vector v = face->ProjectOnBasis(faceBasis, func);

				vectorOfDoFs.segment(i, HHO->nFaceUnknowns) = face->InvFaceMassMatrix()*face->ProjectOnBasis(faceBasis, func);
			}
		);
		return vectorOfDoFs;
	}

	double IntegralFromFaceDoFs(Vector dofs, int polyDegree)
	{
		auto faceBasis = HHO->FaceBasis;

		struct ChunkResult
		{
			double integral = 0;
		};

		ParallelLoop<Face<Dim>*, ChunkResult> parallelLoop(this->_mesh->Faces);
		parallelLoop.Execute([this, dofs, faceBasis, &polyDegree](Face<Dim>* f, ParallelChunk<ChunkResult>* chunk)
			{
				if (f->HasDirichletBC())
					return;
				auto approximate = faceBasis->GetApproximateFunction(dofs, f->Number * faceBasis->Size());
				chunk->Results.integral += f->Integral(approximate, polyDegree);
			});

		double integral = 0;

		parallelLoop.AggregateChunkResults([&integral](ChunkResult chunkResult)
			{
				integral += chunkResult.integral;
			});
		return integral;
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

