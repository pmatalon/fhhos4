#pragma once
#include "../Problem/DiffusionProblem.h"
#include "../Geometry/3D/Tetrahedron.h"
#include "Diff_HHOElement.h"
#include "../Utils/ElementParallelLoop.h"
using namespace std;

template <int Dim>
class Diffusion_HHO : public DiffusionProblem<Dim>
{
private:
	bool _staticCondensation = false;

	// We save what we need to reconstruct the higher-order approximation after solving the linear system
	SparseMatrix Atf;
	Vector bt;
	Vector _dirichletCond;
public:
	HHOParameters<Dim>* HHO;
	Vector ReconstructedSolution;
	Vector GlobalHybridSolution;

	Diffusion_HHO(Mesh<Dim>* mesh, TestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool staticCondensation, string outputDirectory)
		: DiffusionProblem<Dim>(mesh, testCase, outputDirectory)
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

	Diffusion_HHO<Dim>* GetProblemOnCoarserMesh()
	{
		HHOParameters<Dim>* coarseHHO = new HHOParameters<Dim>(this->_mesh->CoarseMesh, HHO->Stabilization, HHO->ReconstructionBasis, HHO->CellBasis, HHO->FaceBasis);
		return new Diffusion_HHO<Dim>(this->_mesh->CoarseMesh, this->_testCase, coarseHHO, _staticCondensation, this->_outputDirectory);
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
		if ((  this->_testCase->Code().compare("sine") == 0
			|| this->_testCase->Code().compare("poly") == 0
			|| this->_testCase->Code().compare("zero") == 0) && this->_diffusionField->IsHomogeneous) // solution is H^{r+2} with r in <= k
		{
			double r = k;
			if (k == 0)
				upperBound = constant * pow(h, 2);
			else
				upperBound = constant * pow(h, r + 2);
		}
		else if (this->_testCase->Code().compare("kellogg") == 0) // solution is H^{1+eps}
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

	void PrintDiscretization()
	{
		cout << "Mesh: " << this->_mesh->Description() << endl;
		cout << "    Elements  : " << HHO->nElements << endl;
		cout << "    Faces     : " << HHO->nFaces << " (" << HHO->nInteriorFaces << " interior + " << HHO->nBoundaryFaces << " boundary)" << endl;
		cout << "    h         : " << scientific << this->_mesh->H() << defaultfloat << endl;
		cout << "    Regularity: " << this->_mesh->Regularity() << defaultfloat << endl;
		cout << "Discretization: Hybrid High-Order (k = " << HHO->FaceBasis->GetDegree() << ")" << endl;
		cout << "    Reconstruction basis: " << HHO->ReconstructionBasis->Name() << endl;
		cout << "    Cell basis          : " << HHO->CellBasis->Name() << endl;
		cout << "    Face basis          : " << HHO->FaceBasis->Name() << endl;
		cout << "Cell unknowns : " << HHO->nTotalCellUnknowns << " (" << HHO->CellBasis->Size() << " per cell)" << endl;
		cout << "Face unknowns : " << HHO->nTotalFaceUnknowns << " (" << HHO->FaceBasis->Size() << " per interior face)" << endl;
		cout << "Total unknowns: " << HHO->nTotalHybridUnknowns << endl;
		cout << "System size   : " << (this->_staticCondensation ? HHO->nTotalFaceUnknowns : HHO->nTotalHybridUnknowns) << " (" << (this->_staticCondensation ? "statically condensed" : "no static condensation") << ")" << endl;
	}

	//--------------------------------------------//
	// Performs the assembly of the linear system //
	//--------------------------------------------//

	void Assemble(ActionsArguments actions) override
	{
		Mesh<Dim>* mesh = this->_mesh;
		FunctionalBasis<Dim>* reconstructionBasis = HHO->ReconstructionBasis;
		FunctionalBasis<Dim>* cellBasis = HHO->CellBasis;
		FunctionalBasis<Dim - 1>* faceBasis = HHO->FaceBasis;

		if (actions.LogAssembly)
		{
			this->PrintPhysicalProblem();
			this->PrintDiscretization();
		}

		string matrixFilePath				= this->GetFilePath("A");
		string consistencyFilePath			= this->GetFilePath("A_cons");
		string stabilizationFilePath		= this->GetFilePath("A_stab");
		string reconstructionMatrixFilePath	= this->GetFilePath("Reconstruct");
		string rhsFilePath					= this->GetFilePath("b");

		if (actions.LogAssembly)
			cout << endl << "Assembly..." << endl;

		Vector globalRHS;
		if (actions.AssembleRightHandSide)
		{
			globalRHS = Vector(HHO->nTotalHybridUnknowns);
			globalRHS.tail(HHO->nTotalFaceUnknowns) = Vector::Zero(HHO->nTotalFaceUnknowns);
		}

		// Compute some useful integrals on reference element and store them
		// - Cartesian element
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreCellMassMatrix(cellBasis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreReconstructMassMatrix(reconstructionBasis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreCellStiffnessMatrix(cellBasis);
		if (this->_diffusionField->K1)
			CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreReconstructK1StiffnessMatrix(this->_diffusionField->K1, reconstructionBasis);
		if (this->_diffusionField->K2)
			CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreReconstructK2StiffnessMatrix(this->_diffusionField->K2, reconstructionBasis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix(cellBasis, reconstructionBasis);
		CartesianShape<Dim, Dim - 1>::InitReferenceShape()->ComputeAndStoreFaceMassMatrix(faceBasis);
		if (Dim == 2)
		{
			// - TriangularElement
			Triangle::InitReferenceShape()->ComputeAndStoreCellMassMatrix((FunctionalBasis<2>*)cellBasis);
			Triangle::InitReferenceShape()->ComputeAndStoreReconstructMassMatrix((FunctionalBasis<2>*)reconstructionBasis);
			//Triangle::InitReferenceShape().ComputeAndStoreCellStiffnessMatrix((FunctionalBasis<2>*)cellBasis);
			//Triangle::InitReferenceShape().ComputeAndStoreReconstructK1StiffnessMatrix(this->_diffusionField->K1, reconstructionBasis);
			//Triangle::InitReferenceShape().ComputeAndStoreReconstructK2StiffnessMatrix(this->_diffusionField->K2, reconstructionBasis);
			Triangle::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix((FunctionalBasis<2>*)cellBasis, (FunctionalBasis<2>*)reconstructionBasis);
		}
		else if (Dim == 3)
		{
			// - TetrahedralElement
			Tetrahedron::InitReferenceShape()->ComputeAndStoreCellMassMatrix((FunctionalBasis<3>*)cellBasis);
			Tetrahedron::InitReferenceShape()->ComputeAndStoreReconstructMassMatrix((FunctionalBasis<3>*)reconstructionBasis);
			Tetrahedron::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix((FunctionalBasis<3>*)cellBasis, (FunctionalBasis<3>*)reconstructionBasis);
		}

		this->InitHHO();

		//-------------------------------//
		// Parallel loop on the elements //
		//-------------------------------//

		struct AssemblyResult
		{
			NonZeroCoefficients MatrixCoeffs;
			NonZeroCoefficients ConsistencyCoeffs;
			NonZeroCoefficients StabilizationCoeffs;
			NonZeroCoefficients ReconstructionCoeffs;
		};

		ParallelLoop<Element<Dim>*, AssemblyResult> parallelLoop(mesh->Elements);

		parallelLoop.InitChunks([this, actions](ParallelChunk<AssemblyResult>* chunk)
			{
				BigNumber nnzApproximate = chunk->Size() * (pow(HHO->nCellUnknowns, 2) + 4 * pow(HHO->nFaceUnknowns, 2) + HHO->nCellUnknowns * 4 * HHO->nFaceUnknowns);
				chunk->Results.MatrixCoeffs         = NonZeroCoefficients(nnzApproximate);
				chunk->Results.ConsistencyCoeffs    = NonZeroCoefficients(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
				chunk->Results.StabilizationCoeffs  = NonZeroCoefficients(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
				chunk->Results.ReconstructionCoeffs = NonZeroCoefficients(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
			});

		parallelLoop.Execute([this, mesh, cellBasis, faceBasis, reconstructionBasis, actions, &globalRHS](Element<Dim>* e, ParallelChunk<AssemblyResult>* chunk)
			{
				Diff_HHOElement<Dim>* element = dynamic_cast<Diff_HHOElement<Dim>*>(e);

				if (!this->_staticCondensation)
				{
					//-----------------------//
					// Att (cell/cell terms) //
					//-----------------------//

					for (BasisFunction<Dim>* cellPhi1 : cellBasis->LocalFunctions)
					{
						BigNumber i = DOFNumber(element, cellPhi1);
						for (BasisFunction<Dim>* cellPhi2 : cellBasis->LocalFunctions)
						{
							BigNumber j = DOFNumber(element, cellPhi2);

							double matrixTerm = element->MatrixTerm(cellPhi1, cellPhi2);
							chunk->Results.MatrixCoeffs.Add(i, j, matrixTerm);
							if (actions.ExportAssemblyTermMatrices)
							{
								double consistencyTerm = element->ConsistencyTerm(cellPhi1, cellPhi2);
								chunk->Results.ConsistencyCoeffs.Add(i, j, consistencyTerm);

								double stabilizationTerm = element->StabilizationTerm(cellPhi1, cellPhi2);
								chunk->Results.StabilizationCoeffs.Add(i, j, stabilizationTerm);
							}
						}
					}
				}

				//----------------------------//
				// Atf, Aft (cell/face terms) //
				//----------------------------//

				for (BasisFunction<Dim>* cellPhi : cellBasis->LocalFunctions)
				{
					BigNumber i = DOFNumber(element, cellPhi);
					for (auto face : element->Faces)
					{
						for (BasisFunction<Dim - 1>* facePhi : faceBasis->LocalFunctions)
						{
							BigNumber j = DOFNumber(face, facePhi);

							double matrixTerm = element->MatrixTerm(face, cellPhi, facePhi);
							chunk->Results.MatrixCoeffs.Add(i, j, matrixTerm);
							chunk->Results.MatrixCoeffs.Add(j, i, matrixTerm);
							if (actions.ExportAssemblyTermMatrices && !face->HasDirichletBC())
							{
								double consistencyTerm = element->ConsistencyTerm(face, cellPhi, facePhi);
								chunk->Results.ConsistencyCoeffs.Add(i, j, consistencyTerm);
								chunk->Results.ConsistencyCoeffs.Add(j, i, consistencyTerm);

								double stabilizationTerm = element->StabilizationTerm(face, cellPhi, facePhi);
								chunk->Results.StabilizationCoeffs.Add(i, j, stabilizationTerm);
								chunk->Results.StabilizationCoeffs.Add(j, i, stabilizationTerm);
							}
						}
					}
				}

				//-----------------------//
				// Aff (face/face terms) //
				//-----------------------//

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
								chunk->Results.MatrixCoeffs.Add(i, j, matrixTerm);
								if (actions.ExportAssemblyTermMatrices && !face1->HasDirichletBC() && !face2->HasDirichletBC())
								{
									double consistencyTerm = element->ConsistencyTerm(face1, facePhi1, face2, facePhi2);
									chunk->Results.ConsistencyCoeffs.Add(i, j, consistencyTerm);

									double stabilizationTerm = element->StabilizationTerm(face1, facePhi1, face2, facePhi2);
									chunk->Results.StabilizationCoeffs.Add(i, j, stabilizationTerm);
								}
							}
						}
					}
				}

				//-----------------//
				// Right-hand side //
				//-----------------//

				if (actions.AssembleRightHandSide)
				{
					for (BasisFunction<Dim>* cellPhi : cellBasis->LocalFunctions)
					{
						BigNumber i = DOFNumber(element, cellPhi);
						globalRHS(i) = element->SourceTerm(cellPhi, this->_sourceFunction);
					}
				}

				//------------------------------------------------//
				// Global reconstruction matrix (only for export) //
				//------------------------------------------------//

				if (actions.ExportAssemblyTermMatrices)
				{
					for (BasisFunction<Dim>* reconstructPhi : reconstructionBasis->LocalFunctions)
					{
						BigNumber i = element->Number * reconstructionBasis->Size() + reconstructPhi->LocalNumber;
						for (BasisFunction<Dim>* cellPhi : cellBasis->LocalFunctions)
						{
							BigNumber j = DOFNumber(element, cellPhi);
							chunk->Results.ReconstructionCoeffs.Add(i, j, element->ReconstructionTerm(reconstructPhi, cellPhi));
						}
						for (auto face : element->Faces)
						{
							for (BasisFunction<Dim - 1>* facePhi : faceBasis->LocalFunctions)
							{
								BigNumber j = DOFNumber(face, facePhi);
								chunk->Results.ReconstructionCoeffs.Add(i, j, element->ReconstructionTerm(reconstructPhi, face, facePhi));
							}
						}
					}
				}
			});

		//------------------------------------//
		// Aggregation of the parallel chunks //
		//------------------------------------//

		cout << "Memory allocation for the global list of coefficients" << endl;

		BigNumber nnzApproximate = mesh->Elements.size() * (pow(HHO->nCellUnknowns, 2) + 4 * pow(HHO->nFaceUnknowns, 2) + HHO->nCellUnknowns * 4 * HHO->nFaceUnknowns);
		NonZeroCoefficients matrixCoeffs(nnzApproximate);
		NonZeroCoefficients consistencyCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
		NonZeroCoefficients stabilizationCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
		NonZeroCoefficients reconstructionCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);

		cout << "Aggregation of the parallel chunks into the global list" << endl;

		parallelLoop.AggregateChunkResults([&matrixCoeffs, &consistencyCoeffs, &stabilizationCoeffs, &reconstructionCoeffs, actions](AssemblyResult& chunkResult)
			{
				matrixCoeffs.Add(chunkResult.MatrixCoeffs);
				chunkResult.MatrixCoeffs.Clear();
				if (actions.ExportAssemblyTermMatrices)
				{
					consistencyCoeffs.Add(chunkResult.ConsistencyCoeffs);
					stabilizationCoeffs.Add(chunkResult.StabilizationCoeffs);
					reconstructionCoeffs.Add(chunkResult.ReconstructionCoeffs);
				}
			});

		//-----------------------------------//
		// Delete now useless local matrices //
		//-----------------------------------//

		cout << "Delete now useless local matrices" << endl;

		ElementParallelLoop<Dim> parallelLoopE(mesh->Elements);
		parallelLoopE.Execute([](Element<Dim>* element)
			{
				Diff_HHOElement<Dim>* e = dynamic_cast<Diff_HHOElement<Dim>*>(element);
				e->DeleteUselessMatricesAfterAssembly();
			});

		FaceParallelLoop<Dim> parallelLoopF(mesh->Faces);
		parallelLoopF.Execute([](Face<Dim>* face)
			{
				Diff_HHOFace<Dim>* f = dynamic_cast<Diff_HHOFace<Dim>*>(face);
				f->DeleteUselessMatricesAfterAssembly();
			});


		cout << "Reserve memory for extendedMatrix" << endl;

		SparseMatrix extendedMatrix = SparseMatrix(HHO->nTotalHybridCoeffs, HHO->nTotalHybridCoeffs);
		extendedMatrix.reserve(matrixCoeffs.Size());

		cout << "Fill extendedMatrix" << endl;
		matrixCoeffs.Fill(extendedMatrix);

		cout << "Clear coefficients" << endl;
		matrixCoeffs.Clear();
		consistencyCoeffs.Clear();
		stabilizationCoeffs.Clear();
		reconstructionCoeffs.Clear();

		// The global system matrix is the extended matrix where the Dirichlet "unknowns" are eliminated
		SparseMatrix globalSystemMatrix = extendedMatrix.topLeftCorner(HHO->nTotalHybridUnknowns, HHO->nTotalHybridUnknowns);

		this->Atf = globalSystemMatrix.topRightCorner(HHO->nTotalCellUnknowns, HHO->nTotalFaceUnknowns); // save this part for the reconstruction

		//---------------------------------//
		// Boundary conditions enforcement //
		//---------------------------------//

		// Neumann
		if (actions.AssembleRightHandSide)
		{
			ParallelLoop<Face<Dim>*>::Execute(this->_mesh->NeumannFaces, [this, faceBasis, &globalRHS](Face<Dim>* f)
				{
					Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);
					BigNumber i = FirstDOFGlobalNumber(f);
					globalRHS.segment(i, HHO->nFaceUnknowns) = face->ProjectOnBasis(faceBasis, this->_boundaryConditions->NeumannFunction);
				}
			);
		}

		// Dirichlet
		this->_dirichletCond = Vector(HHO->nDirichletUnknowns);
		assert(!this->_mesh->DirichletFaces.empty());
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->DirichletFaces, [this, faceBasis](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);
				BigNumber i = FirstDOFGlobalNumber(f) - HHO->nTotalHybridUnknowns;
				this->_dirichletCond.segment(i, HHO->nFaceUnknowns) = face->InvFaceMassMatrix()*face->ProjectOnBasis(faceBasis, this->_boundaryConditions->DirichletFunction);
			}
		);
		if (actions.AssembleRightHandSide)
			globalRHS -= extendedMatrix.topRightCorner(HHO->nTotalHybridUnknowns, HHO->nDirichletUnknowns) * this->_dirichletCond;

		//---------------------//
		// Static condensation //
		//---------------------//

		if (this->_staticCondensation)
		{
			if (actions.LogAssembly)
				cout << "Static condensation..." << endl;

			SparseMatrix Aff = globalSystemMatrix.bottomRightCorner(HHO->nTotalFaceUnknowns, HHO->nTotalFaceUnknowns);

			SparseMatrix inverseAtt = GetInverseAtt();

			// Schur complement
			this->A = (Aff - Atf.transpose() * inverseAtt * Atf).pruned();

			if (actions.AssembleRightHandSide)
			{
				this->bt = globalRHS.head(HHO->nTotalCellUnknowns); // save for reconstruction
				Vector bf = globalRHS.tail(HHO->nTotalFaceUnknowns); // not used in the reconstruction, no need to save it

				// Right-hand side of the condensed system
				this->b = bf - Atf.transpose() * inverseAtt * bt;
			}
		}
		else
		{
			this->A = globalSystemMatrix;
			this->b = globalRHS;
		}

		if (actions.LogAssembly)
			cout << Utils::MatrixInfo(this->A, "A") << endl;

		if (!actions.AssembleRightHandSide)
		{
			// We won't be reconstructing the solution, so ne need to keep that
			Utils::Empty(this->Atf);
			Utils::Empty(this->_dirichletCond);
		}

		//------------------//
		//      Export      //
		//------------------//

		if (actions.ExportLinearSystem)
		{
			cout << "Export:" << endl;
			Eigen::saveMarket(this->A, matrixFilePath);
			cout << "Matrix exported to                " << matrixFilePath << endl;

			Eigen::saveMarketVector(this->b, rhsFilePath);
			cout << "RHS exported to                   " << rhsFilePath << endl;
		}

		if (actions.ExportAssemblyTermMatrices)
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
				Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);
				face->InitHHO(HHO);
			}
		);

		// Init Elements //
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this](Element<Dim>* e)
			{
				Diff_HHOElement<Dim>* element = dynamic_cast<Diff_HHOElement<Dim>*>(e);
				element->InitHHO(HHO);
			}
		);
	}

	//-------------------------------------------------------------------------------------------------//
	// After solving the faces, construction of the higher-order approximation using the reconstructor //
	//-------------------------------------------------------------------------------------------------//
	void ReconstructHigherOrderApproximation()
	{
		Vector globalReconstructedSolution(HHO->nElements * HHO->nReconstructUnknowns);
		this->GlobalHybridSolution = Vector(HHO->nTotalHybridCoeffs);

		if (this->_staticCondensation)
		{
			cout << "Solving cell unknowns..." << endl;
			Vector facesSolution = this->SystemSolution;

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
				Diff_HHOElement<Dim>* element = dynamic_cast<Diff_HHOElement<Dim>*>(e);

				Vector localHybridSolution(HHO->nCellUnknowns + HHO->nFaceUnknowns * element->Faces.size());
				localHybridSolution.head(HHO->nCellUnknowns) = this->GlobalHybridSolution.segment(FirstDOFGlobalNumber(element), HHO->nCellUnknowns);
				for (auto face : element->Faces)
					localHybridSolution.segment(element->FirstDOFNumber(face), HHO->nFaceUnknowns) = this->GlobalHybridSolution.segment(FirstDOFGlobalNumber(face), HHO->nFaceUnknowns);

				Vector localReconstructedSolution = element->Reconstruct(localHybridSolution);
				globalReconstructedSolution.segment(element->Number * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns) = localReconstructedSolution;
			});

		this->ReconstructedSolution = globalReconstructedSolution;
	}

	//---------------------------------------//
	//                Exports                //
	//---------------------------------------//
	void ExtractTraceSystemSolution()
	{
		if (this->_staticCondensation)
			Problem<Dim>::ExportSolutionVector(this->SystemSolution, "Faces");
	}

	void ExtractHybridSolution()
	{
		Problem<Dim>::ExportSolutionVector(this->GlobalHybridSolution, "Hybrid");
	}

	void ExportSolutionVector() override
	{
		Problem<Dim>::ExportSolutionVector(this->ReconstructedSolution);
	}

	void ExportSolutionToGMSH() override
	{
		this->_mesh->ExportSolutionToGMSH(this->HHO->ReconstructionBasis, this->ReconstructedSolution, this->GetFilePathPrefix());
	}

	SparseMatrix GetInverseAtt()
	{
		ElementParallelLoop<Dim> parallelLoop(this->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nCellUnknowns * HHO->nCellUnknowns);
		parallelLoop.Execute([this](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* element = dynamic_cast<Diff_HHOElement<Dim>*>(e);
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
				Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);
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

