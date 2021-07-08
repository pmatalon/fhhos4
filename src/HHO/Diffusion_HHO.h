#pragma once
#include "../Problem/DiffusionProblem.h"
#include "../Geometry/3D/Tetrahedron.h"
#include "Diff_HHOElement.h"
#include "../Utils/ElementParallelLoop.h"
#include "Diffusion_HHOMatrix.h"
using namespace std;

template <int Dim>
class Diffusion_HHO : public DiffusionProblem<Dim>
{
private:
	bool _staticCondensation = true;
	bool _saveMatrixBlocks = false;

	vector<Diff_HHOElement<Dim>> _hhoElements;
	vector<Diff_HHOFace<Dim>> _hhoFaces;
public:
	// Matrix parts
	SparseMatrix A_T_T;
	SparseMatrix A_T_ndF; // used to reconstruct the the higher-order approximation after solving the linear system
	SparseMatrix A_ndF_ndF;
private:
	// Cell part of the right-hand side, used to reconstruct the the higher-order approximation after solving the linear system
	Vector B_T; 
	// Solution on the Dirichlet faces
	Vector x_dF;
public:
	HHOParameters<Dim>* HHO;
	Vector ReconstructedSolution;
	Vector GlobalHybridSolution;

	Diffusion_HHO(Mesh<Dim>* mesh, TestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool staticCondensation, bool saveMatrixBlocks, string outputDirectory)
		: DiffusionProblem<Dim>(mesh, testCase, outputDirectory)
	{	
		this->HHO = hho;
		this->_staticCondensation = staticCondensation;
		this->_saveMatrixBlocks = saveMatrixBlocks;

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
		HHOParameters<Dim>* coarseHHO = new HHOParameters<Dim>(this->_mesh->CoarseMesh, HHO->Stabilization, HHO->ReconstructionBasis, HHO->CellBasis, HHO->FaceBasis, HHO->OrthonormalizeBases);
		return new Diffusion_HHO<Dim>(this->_mesh->CoarseMesh, this->_testCase, coarseHHO, _staticCondensation, _saveMatrixBlocks, this->_outputDirectory);
	}

	Diffusion_HHO<Dim>* GetProblemForLowerDegree()
	{
		// Lower the degree of each basis
		FunctionalBasis<Dim>* reconstructionBasis = HHO->CellBasis;
		FunctionalBasis<Dim>* cellBasis = new FunctionalBasis<Dim>(HHO->CellBasis->CreateSameBasisForLowerDegree());
		FunctionalBasis<Dim-1>* faceBasis = new FunctionalBasis<Dim-1>(HHO->FaceBasis->CreateSameBasisForLowerDegree());

		HHOParameters<Dim>* lowerDegreeHHO = new HHOParameters<Dim>(this->_mesh, HHO->Stabilization, reconstructionBasis, cellBasis, faceBasis, HHO->OrthonormalizeBases);
		return new Diffusion_HHO<Dim>(this->_mesh, this->_testCase, lowerDegreeHHO, _staticCondensation, _saveMatrixBlocks, this->_outputDirectory);
	}

	Diffusion_HHO<Dim>* GetProblemOnCoarserMeshAndLowerDegree()
	{
		// Lower the degree of each basis
		FunctionalBasis<Dim>* reconstructionBasis = HHO->CellBasis;
		FunctionalBasis<Dim>* cellBasis = new FunctionalBasis<Dim>(HHO->CellBasis->CreateSameBasisForLowerDegree());
		FunctionalBasis<Dim - 1>* faceBasis = new FunctionalBasis<Dim - 1>(HHO->FaceBasis->CreateSameBasisForLowerDegree());

		HHOParameters<Dim>* lowerDegreeCoarseMeshHHO = new HHOParameters<Dim>(this->_mesh->CoarseMesh, HHO->Stabilization, reconstructionBasis, cellBasis, faceBasis, HHO->OrthonormalizeBases);
		return new Diffusion_HHO<Dim>(this->_mesh->CoarseMesh, this->_testCase, lowerDegreeCoarseMeshHHO, _staticCondensation, _saveMatrixBlocks, this->_outputDirectory);
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
		cout << "    Face basis          : " << (HHO->OrthonormalizeBases ? "orthonormalized_" : "") << HHO->FaceBasis->Name() << endl;
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

		//------------------------//
		//      Memory usage      //
		//------------------------//

		// Element local matrices
		int nFacesPerElement = Dim == 3 ? 4 : 3;
		int nCoeffs_A_T_T = pow(HHO->nCellUnknowns, 2);
		int nCoeffs_A_T_F = HHO->nCellUnknowns * nFacesPerElement * HHO->nFaceUnknowns;
		int nCoeffs_A_F_F = nFacesPerElement * pow(HHO->nFaceUnknowns, 2);

		int nCoeffs_A_T = nCoeffs_A_T_T + 2 * nCoeffs_A_T_F + nCoeffs_A_F_F;

		int nCoeffs_P = HHO->nReconstructUnknowns * (HHO->nCellUnknowns + nFacesPerElement * HHO->nFaceUnknowns);

		int nCoeffs_localElemMatrices = nCoeffs_A_T; // A
		nCoeffs_localElemMatrices += nCoeffs_A_T; // invA
		nCoeffs_localElemMatrices += nCoeffs_P; // P (reconstructor)
		size_t localElemMatricesMemory = nCoeffs_localElemMatrices * sizeof(double);

		size_t totalElemLocalMatricesMemory = mesh->Elements.size() * localElemMatricesMemory;

		// Face local matrices
		int nCoeffs_Mass = pow(HHO->nFaceUnknowns, 2);
		int nCoeffs_invMass = nCoeffs_Mass;

		int nCoeffs_localFaceMatrices = nCoeffs_Mass + nCoeffs_invMass;
		size_t localFaceMatricesMemory = nCoeffs_localFaceMatrices * sizeof(double);

		size_t totalFaceLocalMatricesMemory = mesh->Faces.size() * localFaceMatricesMemory;


		size_t localNonZeroVectorMemory = nCoeffs_A_T * NonZeroCoefficients::SizeOfNonZero();
		size_t globalNonZeroVectorMemory = mesh->Elements.size() * localNonZeroVectorMemory;

		size_t sparseMatrixMemory = Utils::SparseMatrixMemoryUsage(mesh->Elements.size() * nCoeffs_A_T);

		size_t requiredMemory = totalElemLocalMatricesMemory + totalFaceLocalMatricesMemory + Utils::VectorMemoryUsage(HHO->nTotalHybridUnknowns) + globalNonZeroVectorMemory + sparseMatrixMemory;

		if (actions.LogAssembly)
			cout << "\tRequired memory > " << Utils::MemoryString(requiredMemory) << endl;

		//-----------------------------------------//
		// Compute integrals on reference elements //
		//-----------------------------------------//

		if (actions.InitReferenceShapes)
			InitReferenceShapes();

		//----------------------------//
		//   Compute local matrices   //
		//----------------------------//

		if (actions.LogAssembly)
			cout << "\tCompute local matrices (allocation of " + Utils::MemoryString(totalElemLocalMatricesMemory) + ")" << endl;
		this->InitHHO();

		if (!actions.ExportAssemblyTermMatrices)
		{
			ElementParallelLoop<Dim> parallelLoopE(mesh->Elements);
			parallelLoopE.Execute([this](Element<Dim>* element)
				{
					Diff_HHOElement<Dim>* e = HHOElement(element);
					e->DeleteUselessMatricesAfterAssembly(); // Acons, Astab
				});
		}

		// Notations:
		//
		//  [  A_T_T  |          A_T_F          ] [ x_T ]       [ B_T ]
		//  [ --------|------------------------ ] [-----]       [-----]
		//  [         |                         ] [     ]    =  [     ]
		//  [  <sym>  |          A_F_F          ] [ x_F ]       [  0  ]
		//  [         |                         ] [     ]       [     ]
		//
		//                       <=>             (the faces are split into non-Dirichlet faces (ndF) and Dirichlet faces (dF))
		//
		//  [  A_T_T  |  A_T_ndF   |  A_T_dF    ] [ x_T   ]     [ B_T ]
		//  [ --------|------------|----------- ] [-------]     [-----]
		//  [  <sym>  |  A_ndF_ndF |  A_ndF_dF  ] [ x_ndF ]  =  [  0  ]
		//  [ --------|------------|----------- ] [-------]     [-----]
		//  [  <sym>  |   <sym>    |  A_dF_dF   ] [ x_dF  ]     [  0  ]
		//
		// As x_dF is known, it can be passed to the right-hand side (method by elimination of the Dirichlet conditions):
		//
		//               [  A_T_T  |  A_T_ndF   ] [ x_T   ]     [ B_T        ]          [ A_T_dF   * x_dF ]
		//               [ --------|----------- ] [-------]  =  [------------]     -    [-----------------]
		//               [   sym   |  A_ndF_ndF ] [ x_ndF ]     [ 0; B_neumF ]          [ A_ndF_dF * x_dF ]
		//                                                         ^
		//                                                         |
		//                                                      (interior faces; Neumann faces)


		// Allocation of the cell part of the right-hand side
		if (actions.AssembleRightHandSide)
			this->B_T = Vector(HHO->nTotalCellUnknowns);

		//-------------------------------//
		// Parallel loop on the elements //
		//-------------------------------//

		struct AssemblyResult
		{
			A_T_T_Block<Dim> A_T_T_Coeffs;
			A_T_F_Block<Dim> A_T_F_Coeffs;
			A_F_F_Block<Dim> A_F_F_Coeffs;

			NonZeroCoefficients ConsistencyCoeffs;
			NonZeroCoefficients StabilizationCoeffs;
			NonZeroCoefficients ReconstructionCoeffs;
		};

		ParallelLoop<Element<Dim>*, AssemblyResult> parallelLoop(mesh->Elements);

		if (actions.LogAssembly)
			cout << "\tParallel recovery of the non-zero coefficients in chunks (allocation of " + to_string(parallelLoop.NThreads) + "x" + Utils::MemoryString((mesh->Elements.size() / parallelLoop.NThreads) * localNonZeroVectorMemory) + " = " + Utils::MemoryString(globalNonZeroVectorMemory) + ")" << endl;

		parallelLoop.InitChunks([this, nCoeffs_A_T, nCoeffs_A_T_T, nCoeffs_A_T_F, nCoeffs_A_F_F, actions](ParallelChunk<AssemblyResult>* chunk)
			{
				if (actions.ExportAssemblyTermMatrices)
				{
					BigNumber nnzApproximate = chunk->Size() * nCoeffs_A_T;
					chunk->Results.ConsistencyCoeffs = NonZeroCoefficients(nnzApproximate);
					chunk->Results.StabilizationCoeffs = NonZeroCoefficients(nnzApproximate);
					chunk->Results.ReconstructionCoeffs = NonZeroCoefficients(nnzApproximate);
				}
				if (!this->_staticCondensation || this->_saveMatrixBlocks)
				{
					chunk->Results.A_T_T_Coeffs = A_T_T_Block<Dim>(HHO->nCellUnknowns);
					chunk->Results.A_T_T_Coeffs.Reserve(chunk->Size() * nCoeffs_A_T_T);
				}
				chunk->Results.A_T_F_Coeffs = A_T_F_Block<Dim>(HHO->nCellUnknowns, HHO->nFaceUnknowns);
				chunk->Results.A_T_F_Coeffs.Reserve(chunk->Size() * nCoeffs_A_T_F);

				chunk->Results.A_F_F_Coeffs = A_F_F_Block<Dim>(HHO->nFaceUnknowns);
				chunk->Results.A_F_F_Coeffs.Reserve(chunk->Size() * nCoeffs_A_F_F);
			});

		parallelLoop.Execute([this, mesh, cellBasis, faceBasis, reconstructionBasis, actions](Element<Dim>* e, ParallelChunk<AssemblyResult>* chunk)
			{
				Diff_HHOElement<Dim>* element = HHOElement(e);

				if (!this->_staticCondensation || this->_saveMatrixBlocks)
				{
					//-------------------------//
					// A_T_T (cell/cell terms) //
					//-------------------------//

					A_T_T_Block<Dim>& A_T_T = chunk->Results.A_T_T_Coeffs;

					BigNumber i = A_T_T.FirstRow(element);
					A_T_T.AddBlock(i, i, element->A, 0, 0, HHO->nCellUnknowns, HHO->nCellUnknowns);
					if (actions.ExportAssemblyTermMatrices)
					{
						chunk->Results.ConsistencyCoeffs.AddBlock(i, i, element->Acons, 0, 0, HHO->nCellUnknowns, HHO->nCellUnknowns);
						chunk->Results.StabilizationCoeffs.AddBlock(i, i, element->Astab, 0, 0, HHO->nCellUnknowns, HHO->nCellUnknowns);
					}
				}

				//-------------------------//
				// A_T_F (cell/face terms) //
				//-------------------------//

				A_T_F_Block<Dim>& A_T_F = chunk->Results.A_T_F_Coeffs;

				BigNumber i = A_T_F.FirstRow(element);

				for (auto face : element->Faces)
				{
					BigNumber j = A_T_F.FirstCol(face);
					A_T_F.AddBlock(i, j, element->A, 0, element->FirstDOFNumber(face), HHO->nCellUnknowns, HHO->nFaceUnknowns);

					if (actions.ExportAssemblyTermMatrices && !face->HasDirichletBC())
					{
						chunk->Results.ConsistencyCoeffs.AddBlock(i, j, element->Acons, 0, element->FirstDOFNumber(face), HHO->nCellUnknowns, HHO->nFaceUnknowns);
						chunk->Results.StabilizationCoeffs.AddBlock(i, j, element->Astab, 0, element->FirstDOFNumber(face), HHO->nCellUnknowns, HHO->nFaceUnknowns);
					}
				}

				//-------------------------//
				// A_F_F (face/face terms) //
				//-------------------------//

				A_F_F_Block<Dim>& A_F_F = chunk->Results.A_F_F_Coeffs;

				for (auto face1 : element->Faces)
				{
					BigNumber i = A_F_F.FirstRow(face1);
					for (auto face2 : element->Faces)
					{
						BigNumber j = A_F_F.FirstCol(face2);
						A_F_F.AddBlock(i, j, element->A, element->FirstDOFNumber(face1), element->FirstDOFNumber(face2), HHO->nFaceUnknowns, HHO->nFaceUnknowns);
						if (actions.ExportAssemblyTermMatrices && !face1->HasDirichletBC() && !face2->HasDirichletBC())
						{
							chunk->Results.ConsistencyCoeffs.AddBlock(i, j, element->Acons, element->FirstDOFNumber(face1), element->FirstDOFNumber(face2), HHO->nFaceUnknowns, HHO->nFaceUnknowns);
							chunk->Results.StabilizationCoeffs.AddBlock(i, j, element->Astab, element->FirstDOFNumber(face1), element->FirstDOFNumber(face2), HHO->nFaceUnknowns, HHO->nFaceUnknowns);
						}
					}
				}

				//--------------------------//
				// Right-hand side part B_T //
				//--------------------------//

				if (actions.AssembleRightHandSide)
				{
					for (BasisFunction<Dim>* cellPhi : cellBasis->LocalFunctions)
					{
						BigNumber i = DOFNumber(element, cellPhi);
						this->B_T(i) = element->SourceTerm(cellPhi, this->_sourceFunction);
					}
				}

				//------------------------------------------------//
				// Global reconstruction matrix (only for export) //
				//------------------------------------------------//

				if (actions.ExportAssemblyTermMatrices)
				{
					BigNumber i = element->Number() * reconstructionBasis->Size();
					BigNumber j = FirstDOFGlobalNumber(element);
					chunk->Results.ReconstructionCoeffs.AddBlock(i, j, element->P, 0, 0, reconstructionBasis->Size(), HHO->nCellUnknowns);
					for (auto face : element->Faces)
					{
						j = FirstDOFGlobalNumber(face);
						chunk->Results.ReconstructionCoeffs.AddBlock(i, j, element->P, 0, element->FirstDOFNumber(face), reconstructionBasis->Size(), HHO->nFaceUnknowns);
					}
				}
			});
		
		//-----------------------------------//
		// Delete now useless local matrices //
		//-----------------------------------//

		if (actions.LogAssembly)
			cout << "\tDelete now useless local matrices" << endl;

		ElementParallelLoop<Dim> parallelLoopE(mesh->Elements);
		parallelLoopE.Execute([this](Element<Dim>* element)
			{
				Diff_HHOElement<Dim>* e = HHOElement(element);
				e->DeleteUselessMatricesAfterAssembly();
			});

		FaceParallelLoop<Dim> parallelLoopF(mesh->Faces);
		parallelLoopF.Execute([this](Face<Dim>* face)
			{
				Diff_HHOFace<Dim>* f = HHOFace(face);
				f->DeleteUselessMatricesAfterAssembly();
			});

		//------------------------------------//
		// Aggregation of the parallel chunks //
		//------------------------------------//

		if (actions.LogAssembly)
			cout << "\tAggregation of the parallel chunks into the global list of non-zeroes (allocation of " + Utils::MemoryString(globalNonZeroVectorMemory) + ")" << endl;

		A_T_T_Block<Dim> A_T_T_Coeffs(HHO->nCellUnknowns);
		A_T_F_Block<Dim> A_T_F_Coeffs(HHO->nCellUnknowns, HHO->nFaceUnknowns);
		A_F_F_Block<Dim> A_F_F_Coeffs(HHO->nFaceUnknowns);
		if (!this->_staticCondensation || this->_saveMatrixBlocks)
			A_T_T_Coeffs.Reserve(mesh->Elements.size() * nCoeffs_A_T_T);
		A_T_F_Coeffs.Reserve(mesh->Elements.size() * nCoeffs_A_T_F);
		A_F_F_Coeffs.Reserve(mesh->Faces.size() * nCoeffs_A_F_F);

		BigNumber nnzApproximate = mesh->Elements.size() * nCoeffs_A_T;
		NonZeroCoefficients consistencyCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
		NonZeroCoefficients stabilizationCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);
		NonZeroCoefficients reconstructionCoeffs(actions.ExportAssemblyTermMatrices ? nnzApproximate : 0);

		parallelLoop.AggregateChunkResults([this, &A_T_T_Coeffs, &A_T_F_Coeffs, &A_F_F_Coeffs, &consistencyCoeffs, &stabilizationCoeffs, &reconstructionCoeffs, actions](AssemblyResult& chunkResult)
			{
				if (!this->_staticCondensation || this->_saveMatrixBlocks)
				{
					A_T_T_Coeffs.Add(chunkResult.A_T_T_Coeffs);
					chunkResult.A_T_T_Coeffs = A_T_T_Block<Dim>();
				}
				A_T_F_Coeffs.Add(chunkResult.A_T_F_Coeffs);
				chunkResult.A_T_F_Coeffs = A_T_F_Block<Dim>();

				A_F_F_Coeffs.Add(chunkResult.A_F_F_Coeffs);
				chunkResult.A_F_F_Coeffs = A_F_F_Block<Dim>();

				if (actions.ExportAssemblyTermMatrices)
				{
					consistencyCoeffs.Add(chunkResult.ConsistencyCoeffs);
					stabilizationCoeffs.Add(chunkResult.StabilizationCoeffs);
					reconstructionCoeffs.Add(chunkResult.ReconstructionCoeffs);
				}
			});
		//std::this_thread::sleep_for(std::chrono::seconds(5));

		//-------------------------------------//
		//    Assembly of the sparse matrix    //
		//-------------------------------------//
		
		if (!this->_staticCondensation || this->_saveMatrixBlocks)
		{
			this->A_T_T = SparseMatrix(HHO->nTotalCellUnknowns, HHO->nTotalCellUnknowns);
			A_T_T.reserve(A_T_T_Coeffs.Size());
			A_T_T_Coeffs.Fill(A_T_T);
			A_T_T_Coeffs = A_T_T_Block<Dim>();
		}

		if (actions.LogAssembly)
			cout << "\tReserve memory for the block A_T_F (allocation of " + Utils::MemoryString(Utils::SparseMatrixMemoryUsage(A_T_F_Coeffs.Size())) + ")" << endl;
		SparseMatrix A_T_F(HHO->nTotalCellUnknowns, HHO->nTotalFaceCoeffs);
		A_T_F.reserve(A_T_F_Coeffs.Size());
		if (actions.LogAssembly)
			cout << "\tFill A_T_F with non-zero coefficients" << endl;
		A_T_F_Coeffs.Fill(A_T_F);
		if (actions.LogAssembly)
			cout << "\tFree vector of non-zero coefficients (" + Utils::MemoryString(Utils::SparseMatrixMemoryUsage(A_T_F_Coeffs.Size())) + ")" << endl;
		A_T_F_Coeffs = A_T_F_Block<Dim>();

		this->A_T_ndF = A_T_F.topLeftCorner(HHO->nTotalCellUnknowns, HHO->nTotalFaceUnknowns); // save this part for the reconstruction
		
		if (actions.LogAssembly)
			cout << "\tReserve memory for the block A_F_F (allocation of " + Utils::MemoryString(Utils::SparseMatrixMemoryUsage(A_F_F_Coeffs.Size())) + " for " + to_string(A_F_F_Coeffs.Size()) + " non-zeroes)" << endl;
		if (actions.LogAssembly)
			cout << "\tA_F_F(" + to_string(HHO->nTotalFaceCoeffs) + ", " + to_string(HHO->nTotalFaceCoeffs) + ")" << endl;
		SparseMatrix A_F_F(HHO->nTotalFaceCoeffs, HHO->nTotalFaceCoeffs);
		A_F_F.reserve(A_F_F_Coeffs.Size());
		if (actions.LogAssembly)
			cout << "\tFill A_F_F with non-zero coefficients" << endl;
		A_F_F_Coeffs.Fill(A_F_F);
		if (actions.LogAssembly)
			cout << "\tFree vector of non-zero coefficients (" + Utils::MemoryString(Utils::SparseMatrixMemoryUsage(A_F_F_Coeffs.Size())) + ")" << endl;
		A_F_F_Coeffs = A_F_F_Block<Dim>();

		//---------------------------------//
		// Boundary conditions enforcement //
		//---------------------------------//

		// 0 on the interior faces, computation for Neumann faces
		Vector B_ndF = Vector::Zero(HHO->nTotalFaceUnknowns);
		if (actions.AssembleRightHandSide)
		{
			ParallelLoop<Face<Dim>*>::Execute(this->_mesh->NeumannFaces, [this, faceBasis, &B_ndF](Face<Dim>* f)
				{
					Diff_HHOFace<Dim>* face = HHOFace(f);
					BigNumber i = FirstDOFGlobalNumber(face) - HHO->nTotalCellUnknowns;
					B_ndF.segment(i, HHO->nFaceUnknowns) = face->ProjectOnBasis(faceBasis, this->_boundaryConditions->NeumannFunction);
				}
			);
		}

		// Solution on the Dirichlet faces
		this->x_dF = Vector(HHO->nDirichletCoeffs);
		assert(!this->_mesh->DirichletFaces.empty());
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->DirichletFaces, [this, faceBasis](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = FirstDOFGlobalNumber(face) - HHO->nTotalHybridUnknowns;
				this->x_dF.segment(i, HHO->nFaceUnknowns) = face->InvMassMatrix()*face->ProjectOnBasis(faceBasis, this->_boundaryConditions->DirichletFunction);
			}
		);

		// Update the right-hand side
		if (actions.AssembleRightHandSide)
		{
			SparseMatrix A_ndF_dF = A_F_F.topRightCorner(HHO->nTotalFaceUnknowns, HHO->nDirichletCoeffs);
			SparseMatrix A_T_dF   = A_T_F.topRightCorner(HHO->nTotalCellUnknowns, HHO->nDirichletCoeffs);
			Utils::Empty(A_T_F);

			this->B_T   -= A_T_dF   * this->x_dF; // save for reconstruction
			      B_ndF -= A_ndF_dF * this->x_dF; // not used in the reconstruction, no need to save it
		}

		

		//---------------------//
		// Static condensation //
		//---------------------//

		this->A_ndF_ndF = A_F_F.topLeftCorner(HHO->nTotalFaceUnknowns, HHO->nTotalFaceUnknowns);
		Utils::Empty(A_F_F);

		if (this->_staticCondensation)
		{
			if (actions.LogAssembly)
				cout << "Static condensation..." << endl;

			SparseMatrix inv_A_T_T = Inverse_A_T_T();

			// Schur complement
			Problem<Dim>::A = (A_ndF_ndF - A_T_ndF.transpose() * inv_A_T_T * A_T_ndF).pruned();

			if (actions.AssembleRightHandSide)
				// Right-hand side of the condensed system
				Problem<Dim>::b = B_ndF - A_T_ndF.transpose() * inv_A_T_T * B_T;
		}
		else
		{
			Problem<Dim>::A = SparseMatrix(HHO->nTotalHybridUnknowns, HHO->nTotalHybridUnknowns);
			NonZeroCoefficients Acoeffs(A_T_T.nonZeros() + 2 * A_T_ndF.nonZeros() + A_ndF_ndF.nonZeros());
			Acoeffs.Add(           0,            0, A_T_T);    // topLeftCorner
			Acoeffs.Add(           0, A_T_T.cols(), A_T_ndF);  // topRightCorner
			Acoeffs.Add(A_T_T.rows(),            0, A_T_ndF.transpose().eval()); // bottomLeftCorner
			Acoeffs.Add(A_T_T.rows(), A_T_T.cols(), A_ndF_ndF); // bottomRightCorner
			Acoeffs.Fill(Problem<Dim>::A);

			Problem<Dim>::b = Vector(HHO->nTotalHybridUnknowns);
			Problem<Dim>::b.head(B_T.rows()) = B_T;
			Problem<Dim>::b.tail(B_ndF.rows()) = B_ndF;
		}

		if (actions.LogAssembly)
			cout << Utils::MatrixInfo(this->A, "A") << endl;

		if (this->_staticCondensation && !this->_saveMatrixBlocks)
		{
			Utils::Empty(this->A_T_T);
			Utils::Empty(this->A_ndF_ndF);
		}

		if (!actions.AssembleRightHandSide)
		{
			// We won't be reconstructing the solution, so no need to keep those
			Utils::Empty(this->A_T_ndF);
			Utils::Empty(this->B_T);
			Utils::Empty(this->x_dF);
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


	// Compute some useful integrals on reference element and store them
	void InitReferenceShapes()
	{
		FunctionalBasis<Dim>* reconstructionBasis = HHO->ReconstructionBasis;
		FunctionalBasis<Dim>* cellBasis = HHO->CellBasis;
		FunctionalBasis<Dim - 1>* faceBasis = HHO->FaceBasis;

		// - Cartesian element
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreMassMatrix(cellBasis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreMassMatrix(reconstructionBasis);
		if (this->_diffusionField->K1)
			CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreReconstructK1StiffnessMatrix(this->_diffusionField->K1, reconstructionBasis);
		if (this->_diffusionField->K2)
			CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreReconstructK2StiffnessMatrix(this->_diffusionField->K2, reconstructionBasis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix(cellBasis, reconstructionBasis);
		// - Cartesian face
		CartesianShape<Dim, Dim - 1>::InitReferenceShape()->ComputeAndStoreMassMatrix(faceBasis);
		if (Dim == 2)
		{
			// - Triangular element
			Triangle::InitReferenceShape()->ComputeAndStoreMassMatrix((FunctionalBasis<2>*)cellBasis);
			Triangle::InitReferenceShape()->ComputeAndStoreMassMatrix((FunctionalBasis<2>*)reconstructionBasis);
			Triangle::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix((FunctionalBasis<2>*)cellBasis, (FunctionalBasis<2>*)reconstructionBasis);
		}
		else if (Dim == 3)
		{
			// - Tetrahedral element
			Tetrahedron::InitReferenceShape()->ComputeAndStoreMassMatrix((FunctionalBasis<3>*)cellBasis);
			Tetrahedron::InitReferenceShape()->ComputeAndStoreMassMatrix((FunctionalBasis<3>*)reconstructionBasis);
			Tetrahedron::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix((FunctionalBasis<3>*)cellBasis, (FunctionalBasis<3>*)reconstructionBasis);
			// - Triangular face
			Triangle::InitReferenceShape()->ComputeAndStoreMassMatrix((FunctionalBasis<2>*)faceBasis);
		}
	}


	//------------------------------//
	//           Init HHO           //
	//------------------------------//
	void InitHHO()
	{
		InitHHO_Faces();
		InitHHO_Elements();
	}
	void InitHHO_Faces()
	{
		if (!_hhoFaces.empty())
			Utils::Warning("The HHO faces of this problem have already been initialized.");

		_hhoFaces = vector<Diff_HHOFace<Dim>>(this->_mesh->Faces.size());

		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->Faces, [this](Face<Dim>* f)
			{
				_hhoFaces[f->Number].MeshFace = f;
				_hhoFaces[f->Number].InitHHO(HHO);
			}
		);
	}
	void InitHHO_Elements()
	{
		if (!_hhoElements.empty())
			Utils::Warning("The HHO elements of this problem have already been initialized.");

		_hhoElements = vector<Diff_HHOElement<Dim>>(this->_mesh->Elements.size());

		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this](Element<Dim>* e)
			{
				_hhoElements[e->Number].MeshElement = e;

				for (Face<Dim>* f : e->Faces)
					_hhoElements[e->Number].Faces.push_back(&this->_hhoFaces[f->Number]);

				_hhoElements[e->Number].InitHHO(HHO);
			}
		);
	}

	Diff_HHOElement<Dim>* HHOElement(Element<Dim>* e)
	{
		return &this->_hhoElements[e->Number];
	}
	Diff_HHOFace<Dim>* HHOFace(Face<Dim>* f)
	{
		return &this->_hhoFaces[f->Number];
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

			SparseMatrix inverseAtt = Inverse_A_T_T();

			this->GlobalHybridSolution.head(HHO->nTotalCellUnknowns) = inverseAtt * (B_T - A_T_ndF * facesSolution);
			this->GlobalHybridSolution.segment(HHO->nTotalCellUnknowns, HHO->nTotalFaceUnknowns) = facesSolution;
		}
		else
			this->GlobalHybridSolution.head(HHO->nTotalHybridUnknowns) = this->SystemSolution;

		// Dirichlet boundary conditions
		this->GlobalHybridSolution.tail(HHO->nDirichletCoeffs) = this->x_dF;



		cout << "Reconstruction of higher order approximation..." << endl;
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &globalReconstructedSolution](Element<Dim>* e)
			{
				HHOParameters<Dim>* HHO = this->HHO;
				Diff_HHOElement<Dim>* element = HHOElement(e);

				Vector localHybridSolution(HHO->nCellUnknowns + HHO->nFaceUnknowns * element->Faces.size());
				localHybridSolution.head(HHO->nCellUnknowns) = this->GlobalHybridSolution.segment(FirstDOFGlobalNumber(element), HHO->nCellUnknowns);
				for (auto face : element->Faces)
					localHybridSolution.segment(element->FirstDOFNumber(face), HHO->nFaceUnknowns) = this->GlobalHybridSolution.segment(FirstDOFGlobalNumber(face), HHO->nFaceUnknowns);

				Vector localReconstructedSolution = element->Reconstruct(localHybridSolution);
				globalReconstructedSolution.segment(element->Number() * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns) = localReconstructedSolution;
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
		this->_mesh->ExportToGMSH(this->HHO->ReconstructionBasis, this->ReconstructedSolution, this->GetFilePathPrefix(), "potential");
	}

	void ExportErrorToGMSH(const Vector& faceCoeffs) override
	{
		SparseMatrix inverseAtt = Inverse_A_T_T();
		Vector cellCoeffs = inverseAtt * (B_T - A_T_ndF * faceCoeffs);
		this->_mesh->ExportToGMSH(this->HHO->CellBasis, cellCoeffs, this->GetFilePathPrefix(), "error");
	}

	SparseMatrix Inverse_A_T_T()
	{
		ElementParallelLoop<Dim> parallelLoop(this->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nCellUnknowns * HHO->nCellUnknowns);
		parallelLoop.Execute([this](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* element = HHOElement(e);
				chunk->Results.Coeffs.Add(element->Number() * HHO->nCellUnknowns, element->Number() * HHO->nCellUnknowns, element->invAtt);
			});
		SparseMatrix invA_T_T = SparseMatrix(HHO->nTotalCellUnknowns, HHO->nTotalCellUnknowns);
		parallelLoop.Fill(invA_T_T);
		return invA_T_T;
	}

	Vector ProjectOnFaceDiscreteSpace(DomFunction func)
	{
		Vector vectorOfDoFs = Vector(HHO->nTotalFaceUnknowns);
		auto faceBasis = HHO->FaceBasis;
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->Faces, [this, faceBasis, &vectorOfDoFs, func](Face<Dim>* f)
			{
				if (f->HasDirichletBC())
					return;
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = face->Number() * faceBasis->Size();// FirstDOFGlobalNumber(f);

				DenseMatrix m = face->InvMassMatrix();
				Vector v = face->ProjectOnBasis(faceBasis, func);

				vectorOfDoFs.segment(i, HHO->nFaceUnknowns) = face->InvMassMatrix()*face->ProjectOnBasis(faceBasis, func);
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
	BigNumber DOFNumber(Diff_HHOElement<Dim>* element, BasisFunction<Dim>* cellPhi)
	{
		return FirstDOFGlobalNumber(element) + cellPhi->LocalNumber;
	}
	BigNumber FirstDOFGlobalNumber(Diff_HHOElement<Dim>* element)
	{
		return element->Number() * HHO->CellBasis->Size();
	}
	BigNumber FirstDOFGlobalNumber(Diff_HHOFace<Dim>* face)
	{
		return this->HHO->nTotalCellUnknowns + face->Number() * HHO->FaceBasis->Size();
	}
};

