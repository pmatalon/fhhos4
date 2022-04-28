#pragma once
#include "../../Mesh/Mesh.h"
#include "../../Utils/Utils.h"
#ifdef ENABLE_3D
#include "../../Geometry/3D/Tetrahedron.h"
#endif // ENABLE_3D
#include "../../Geometry/CartesianShape.h"
#include "../../Geometry/2D/Triangle.h"
#include "DiscreteSpaces/HHOSkeletonSpace.h"
#include "DiscreteSpaces/HHOReconstructSpace.h"
#include "DiscreteSpaces/HHOBoundarySpace.h"
#include "Diffusion_HHOMatrix.h"
using namespace std;

template <int Dim>
class Diffusion_HHO
{
private:
	string _filePrefix;
	bool _staticCondensation = true;
	bool _saveMatrixBlocks = false;

	vector<Diff_HHOElement<Dim>> _hhoElements;
	vector<Diff_HHOFace<Dim>> _hhoFaces;
public:
	Mesh<Dim>* _mesh;
	DiffusionTestCase<Dim>* TestCase;
	SparseMatrix A;
	Vector b;

	HHOSkeletonSpace<Dim> SkeletonSpace;
	HHOReconstructSpace<Dim> ReconstructSpace;
	HHOBoundarySpace<Dim> BoundarySpace;

	// Matrix parts:
	// Used in reconstruction: A_T_ndF
	// Used in right-hand side: A_T_dF, A_ndF_dF
	SparseMatrix A_T_T, A_T_ndF,   A_T_dF;
	SparseMatrix        A_ndF_ndF, A_ndF_dF;
private:
	// Solution on the Dirichlet faces
	Vector x_dF;
	// Part of the right-hand side corresponding to the cells (must be saved because used to reconstruct the the higher-order approximation after solving the linear system)
	Vector B_T;
	// Part of the right-hand side corresponding to the faces
	Vector B_ndF;
	// Part of the right-hand side corresponding to Neumann faces
	//Vector B_neumF;
public:
	HHOParameters<Dim>* HHO;

	Diffusion_HHO()
	{}

	Diffusion_HHO(Mesh<Dim>* mesh, DiffusionTestCase<Dim>* testCase, HHOParameters<Dim>* hho, bool staticCondensation, bool saveMatrixBlocks)
	{	
		this->_mesh = mesh;
		this->TestCase = testCase;
		this->HHO = hho;
		this->_staticCondensation = staticCondensation;
		this->_saveMatrixBlocks = saveMatrixBlocks;

		//Problem<Dim>::AddFilePrefix("_HHO_" + HHO->ReconstructionBasis->Name() + (_staticCondensation ? "" : "_nostaticcond"));

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
		HHOParameters<Dim>* coarseHHO = new HHOParameters<Dim>(this->_mesh->CoarseMesh, HHO->Stabilization, HHO->ReconstructionBasis, HHO->CellBasis, HHO->FaceBasis, HHO->OrthogonalizeElemBasesCode, HHO->OrthogonalizeFaceBasesCode);
		return new Diffusion_HHO<Dim>(this->_mesh->CoarseMesh, this->TestCase, coarseHHO, _staticCondensation, _saveMatrixBlocks);
	}

	Diffusion_HHO<Dim>* GetProblemForLowerDegree(int faceDegree)
	{
		// Lower the degree of each basis
		int k = HHO->FaceBasis->GetDegree();
		int reduction = k - faceDegree;
		int newReconstructDegree = HHO->ReconstructionBasis->GetDegree() - reduction;
		int newCellDegree = HHO->CellBasis->GetDegree() - reduction;
		FunctionalBasis<Dim>* reconstructionBasis = new FunctionalBasis<Dim>(HHO->ReconstructionBasis->CreateSameBasisForDegree(newReconstructDegree));
		FunctionalBasis<Dim>* cellBasis = new FunctionalBasis<Dim>(HHO->CellBasis->CreateSameBasisForDegree(newCellDegree));
		FunctionalBasis<Dim-1>* faceBasis = new FunctionalBasis<Dim-1>(HHO->FaceBasis->CreateSameBasisForDegree(faceDegree));

		HHOParameters<Dim>* lowerDegreeHHO = new HHOParameters<Dim>(this->_mesh, HHO->Stabilization, reconstructionBasis, cellBasis, faceBasis, HHO->OrthogonalizeElemBasesCode, HHO->OrthogonalizeFaceBasesCode);
		return new Diffusion_HHO<Dim>(this->_mesh, this->TestCase, lowerDegreeHHO, _staticCondensation, _saveMatrixBlocks);
	}

	Diffusion_HHO<Dim>* GetProblemOnCoarserMeshAndLowerDegree(int faceDegree)
	{
		// Lower the degree of each basis
		int k = HHO->FaceBasis->GetDegree();
		int reduction = k - faceDegree;
		int newReconstructDegree = HHO->ReconstructionBasis->GetDegree() - reduction;
		int newCellDegree = HHO->CellBasis->GetDegree() - reduction;
		FunctionalBasis<Dim>* reconstructionBasis = new FunctionalBasis<Dim>(HHO->ReconstructionBasis->CreateSameBasisForDegree(newReconstructDegree));
		FunctionalBasis<Dim>* cellBasis = new FunctionalBasis<Dim>(HHO->CellBasis->CreateSameBasisForDegree(newCellDegree));
		FunctionalBasis<Dim - 1>* faceBasis = new FunctionalBasis<Dim - 1>(HHO->FaceBasis->CreateSameBasisForDegree(faceDegree));

		HHOParameters<Dim>* lowerDegreeCoarseMeshHHO = new HHOParameters<Dim>(this->_mesh->CoarseMesh, HHO->Stabilization, reconstructionBasis, cellBasis, faceBasis, HHO->OrthogonalizeElemBasesCode, HHO->OrthogonalizeFaceBasesCode);
		return new Diffusion_HHO<Dim>(this->_mesh->CoarseMesh, this->TestCase, lowerDegreeCoarseMeshHHO, _staticCondensation, _saveMatrixBlocks);
	}

	double L2Error(DomFunction exactSolution, const Vector& reconstructedSolution)
	{
		struct ChunkResult
		{
			double absoluteError = 0;
			double normExactSolution = 0;
		};

		ParallelLoop<Element<Dim>*, ChunkResult> parallelLoop(this->_mesh->Elements);
		parallelLoop.Execute([this, exactSolution, &reconstructedSolution](Element<Dim>* e, ParallelChunk<ChunkResult>* chunk)
			{
				auto approximate = HHOElement(e)->ReconstructionBasis->GetApproximateFunction(reconstructedSolution, e->Number * HHO->nReconstructUnknowns);
				chunk->Results.absoluteError += e->L2ErrorPow2(approximate, exactSolution);
				chunk->Results.normExactSolution += e->Integral([exactSolution](const DomPoint& p) { return pow(exactSolution(p), 2); });
			});


		double absoluteError = 0;
		double normExactSolution = 0;

		parallelLoop.AggregateChunkResults([&absoluteError, &normExactSolution](ChunkResult& chunkResult)
			{
				absoluteError += chunkResult.absoluteError;
				normExactSolution += chunkResult.normExactSolution;
			});

		absoluteError = sqrt(absoluteError);
		normExactSolution = sqrt(normExactSolution);
		return normExactSolution != 0 ? absoluteError / normExactSolution : absoluteError;
	}

	void AssertSchemeConvergence(double l2Error)
	{
		if (Dim == 1)
			return;

		double h = this->_mesh->H();
		int k = this->HHO->FaceBasis->GetDegree();

		double constant = ConvergenceHiddenConstant(k);
		double upperBound;
		if ((  TestCase->Code().compare("sine") == 0
			|| TestCase->Code().compare("poly") == 0
			|| TestCase->Code().compare("zero") == 0) && TestCase->DiffField.IsHomogeneous) // solution is H^{r+2} with r in <= k
		{
			double r = k;
			if (k == 0)
				upperBound = constant * pow(h, 2);
			else
				upperBound = constant * pow(h, r + 2);
		}
		else if (this->TestCase->Code().compare("kellogg") == 0) // solution is H^{1+eps}
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
		cout << "    h         : " << scientific << this->_mesh->H() << " (average = " << this->_mesh->AverageH() << ")" << defaultfloat << endl;
		cout << "    Regularity: " << this->_mesh->Regularity() << defaultfloat << endl;
		cout << "Discretization: Hybrid High-Order (k = " << HHO->FaceBasis->GetDegree() << ")" << endl;
		cout << "    Reconstruction basis: " << (HHO->OrthogonalizeElemBases() ? (HHO->OrthonormalizeElemBases() ? "orthonormalized_" : "orthogonalized_") : "") << HHO->ReconstructionBasis->Name() << endl;
		cout << "    Cell basis          : " << (HHO->OrthogonalizeElemBases() ? (HHO->OrthonormalizeElemBases() ? "orthonormalized_" : "orthogonalized_") : "") << HHO->CellBasis->Name() << endl;
		cout << "    Face basis          : " << (HHO->OrthogonalizeFaceBases() ? (HHO->OrthonormalizeFaceBases() ? "orthonormalized_" : "orthogonalized_") : "") << HHO->FaceBasis->Name() << endl;
		cout << "Cell unknowns : " << HHO->nTotalCellUnknowns << " (" << HHO->CellBasis->Size() << " per cell)" << endl;
		cout << "Face unknowns : " << HHO->nTotalFaceUnknowns << " (" << HHO->FaceBasis->Size() << " per interior face)" << endl;
		cout << "Total unknowns: " << HHO->nTotalHybridUnknowns << endl;
		cout << "System size   : " << (this->_staticCondensation ? HHO->nTotalFaceUnknowns : HHO->nTotalHybridUnknowns) << " (" << (this->_staticCondensation ? "statically condensed" : "no static condensation") << ")" << endl;
	}

	//--------------------------------------------//
	// Performs the assembly of the linear system //
	//--------------------------------------------//

	void Assemble(const ActionsArguments& actions)
	{
		assert(!actions.Export.LinearSystem && !actions.Export.AssemblyTermMatrices);
		Assemble(actions, ExportModule(""));
	}

	void Assemble(const ActionsArguments& actions, const ExportModule& out)
	{
		Mesh<Dim>* mesh = this->_mesh;

		if (actions.LogAssembly)
			this->PrintDiscretization();

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

		if (!actions.Export.AssemblyTermMatrices)
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
		//  [  A_T_T  |          A_T_F          ] [ x_T ]       [ b_source ]
		//  [ --------|------------------------ ] [-----]       [----------]
		//  [         |                         ] [     ]    =  [          ]
		//  [  <sym>  |          A_F_F          ] [ x_F ]       [ b_F      ]
		//  [         |                         ] [     ]       [          ]
		//
		//                       <=>             (the faces are split into non-Dirichlet faces (ndF) and Dirichlet faces (dF))
		//
		//  [  A_T_T  |  A_T_ndF   |  A_T_dF    ] [ x_T   ]     [ b_source  ]
		//  [ --------|------------|----------- ] [-------]     [-----------]
		//  [  <sym>  |  A_ndF_ndF |  A_ndF_dF  ] [ x_ndF ]  =  [ b_neumann ]
		//  [ --------|------------|----------- ] [-------]     [-----------]
		//  [  <sym>  |   <sym>    |  A_dF_dF   ] [ x_dF  ]     [     0     ]
		// 
		//                                                                   [  0   ]   <---- interior faces
		//                                                 where b_neumann = [------]
		//                                                                   [ neum ]   <---- Neumann faces
		//
		// As x_dF is known, it can be passed to the right-hand side (method by elimination of the Dirichlet conditions):
		//
		//               [  A_T_T  |  A_T_ndF   ] [ x_T   ]     [ b_source  ]       [ A_T_dF   * x_dF ]
		//               [ --------|----------- ] [-------]  =  [-----------]   -   [-----------------]
		//               [   sym   |  A_ndF_ndF ] [ x_ndF ]     [ b_neumann ]       [ A_ndF_dF * x_dF ]
		//                                                     |_____________________________________|
		//                                                                       =
		//                                                                   [ B_T   ]
		//                                                                   [-------]
		//                                                                   [ B_ndF ]
		//
		// Static condensation:
		// 
		//       (A_ndF_ndF - A_T_ndF^T * A_T_T^-1 * A_T_ndF) x_ndF  =  B_ndF - A_T_ndF^T * A_T_T^-1 * B_T
		// 
		// Recover cell unknowns:
		// 
		//        x_T = A_T_T^-1 (B_T - A_T_ndF*x_ndF)


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
				if (actions.Export.AssemblyTermMatrices)
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

		parallelLoop.Execute([this, mesh, actions](Element<Dim>* e, ParallelChunk<AssemblyResult>* chunk)
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
					if (actions.Export.AssemblyTermMatrices)
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

					if (actions.Export.AssemblyTermMatrices && !face->HasDirichletBC())
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
						if (actions.Export.AssemblyTermMatrices && !face1->HasDirichletBC() && !face2->HasDirichletBC())
						{
							chunk->Results.ConsistencyCoeffs.AddBlock(i, j, element->Acons, element->FirstDOFNumber(face1), element->FirstDOFNumber(face2), HHO->nFaceUnknowns, HHO->nFaceUnknowns);
							chunk->Results.StabilizationCoeffs.AddBlock(i, j, element->Astab, element->FirstDOFNumber(face1), element->FirstDOFNumber(face2), HHO->nFaceUnknowns, HHO->nFaceUnknowns);
						}
					}
				}

				//------------------------------------------------//
				// Global reconstruction matrix (only for export) //
				//------------------------------------------------//

				if (actions.Export.AssemblyTermMatrices)
				{
					BigNumber i = element->Number() * HHO->nReconstructUnknowns;
					BigNumber j = FirstDOFGlobalNumber(element);
					chunk->Results.ReconstructionCoeffs.AddBlock(i, j, element->P, 0, 0, HHO->nReconstructUnknowns, HHO->nCellUnknowns);
					for (auto face : element->Faces)
					{
						j = FirstDOFGlobalNumber(face);
						chunk->Results.ReconstructionCoeffs.AddBlock(i, j, element->P, 0, element->FirstDOFNumber(face), HHO->nReconstructUnknowns, HHO->nFaceUnknowns);
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

		if (this->TestCase->BC.Type != PbBoundaryConditions::FullNeumann && Utils::ProgramArgs.Problem.Equation == EquationType::Diffusion)
		{
			// TODO: in the coarse level of multigrid, this should be executed as well
			FaceParallelLoop<Dim> parallelLoopF(mesh->Faces);
			parallelLoopF.Execute([this](Face<Dim>* face)
				{
					Diff_HHOFace<Dim>* f = HHOFace(face);
					f->DeleteUselessMatricesAfterAssembly();
				});
		}

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
		NonZeroCoefficients consistencyCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);
		NonZeroCoefficients stabilizationCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);
		NonZeroCoefficients reconstructionCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);

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

				if (actions.Export.AssemblyTermMatrices)
				{
					consistencyCoeffs.Add(chunkResult.ConsistencyCoeffs);
					stabilizationCoeffs.Add(chunkResult.StabilizationCoeffs);
					reconstructionCoeffs.Add(chunkResult.ReconstructionCoeffs);
				}
			});

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
		this->A_T_dF  = A_T_F.topRightCorner(HHO->nTotalCellUnknowns, HHO->nDirichletCoeffs);
		Utils::Empty(A_T_F);
		
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

		this->A_ndF_dF = A_F_F.topRightCorner(HHO->nTotalFaceUnknowns, HHO->nDirichletCoeffs);

		//-----------------//
		// Right-hand side //
		//-----------------//

		if (actions.AssembleRightHandSide)
		{
			Vector b_source = AssembleSourceTerm(TestCase->SourceFunction);
			this->x_dF  = std::move(AssembleDirichletTerm(TestCase->BC.DirichletFunction));
			Vector b_neumann = AssembleNeumannTerm(TestCase->BC.NeumannFunction);

			// Elimination of the Dirichlet unknowns x_dF from the system
			this->B_T   = ComputeB_T  (b_source,  this->x_dF); // b_source  - A_T_dF   * this->x_dF; // save for reconstruction
			this->B_ndF = ComputeB_ndF(b_neumann, this->x_dF); // b_neumann - A_ndF_dF * this->x_dF;
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

			// Schur complement
			SparseMatrix tmp = SparseMatrix(A_T_ndF.transpose()) * Solve_A_T_T(A_T_ndF);
			this->A = (A_ndF_ndF - tmp).pruned();

			if (actions.AssembleRightHandSide)
				this->b = CondensedRHS(B_T, B_ndF); // B_ndF - A_T_ndF.transpose() * Solve_A_T_T(B_T);
		}
		else
		{
			this->A = SparseMatrix(HHO->nTotalHybridUnknowns, HHO->nTotalHybridUnknowns);
			NonZeroCoefficients Acoeffs(A_T_T.nonZeros() + 2 * A_T_ndF.nonZeros() + A_ndF_ndF.nonZeros());
			Acoeffs.Add(           0,            0, A_T_T);    // topLeftCorner
			Acoeffs.Add(           0, A_T_T.cols(), A_T_ndF);  // topRightCorner
			Acoeffs.Add(A_T_T.rows(),            0, A_T_ndF.transpose().eval()); // bottomLeftCorner
			Acoeffs.Add(A_T_T.rows(), A_T_T.cols(), A_ndF_ndF); // bottomRightCorner
			Acoeffs.Fill(this->A);

			this->b = Vector(HHO->nTotalHybridUnknowns);
			this->b.head(B_T.rows()) = B_T;
			this->b.tail(B_ndF.rows()) = B_ndF;
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
			if (Utils::ProgramArgs.Problem.Equation != EquationType::BiHarmonic)
			{
				// Needed to build the right-hand side separately
				Utils::Empty(this->A_T_dF);
				Utils::Empty(this->A_ndF_dF);
			}
			Utils::Empty(this->B_T);
			Utils::Empty(this->B_ndF);
			Utils::Empty(this->x_dF);
		}

		//------------------//
		//      Export      //
		//------------------//

		if (actions.Export.LinearSystem)
		{
			cout << "Export linear system..." << endl;
			out.ExportMatrix(this->A, "A");
			out.ExportVector(this->b, "b");
		}

		if (actions.Export.AssemblyTermMatrices)
		{
			SparseMatrix Acons = SparseMatrix(HHO->nTotalHybridUnknowns, HHO->nTotalHybridUnknowns);
			consistencyCoeffs.Fill(Acons);

			SparseMatrix Astab = SparseMatrix(HHO->nTotalHybridUnknowns, HHO->nTotalHybridUnknowns);
			stabilizationCoeffs.Fill(Astab);

			SparseMatrix reconstructionMatrix = SparseMatrix(HHO->nElements * HHO->nReconstructUnknowns, HHO->nTotalHybridCoeffs);
			reconstructionCoeffs.Fill(reconstructionMatrix);

			out.ExportMatrix(Acons, "Acons");
			out.ExportMatrix(Astab, "Astab");
			out.ExportMatrix(reconstructionMatrix, "Reconstruct");
		}
	}

	void InitReferenceShapes()
	{
		InitReferenceShapes(this->HHO, &this->TestCase->DiffField);
	}

	// Compute some useful integrals on reference element and store them
	static void InitReferenceShapes(HHOParameters<Dim>* hho, DiffusionField<Dim>* diffField)
	{
		if (hho->OrthogonalizeElemBases() && hho->OrthogonalizeFaceBases())
			return;

		FunctionalBasis<Dim>* reconstructionBasis = hho->ReconstructionBasis;
		FunctionalBasis<Dim>* cellBasis = hho->CellBasis;
		FunctionalBasis<Dim - 1>* faceBasis = hho->FaceBasis;

		// - Cartesian element
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreMassMatrix(cellBasis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreMassMatrix(reconstructionBasis);
		if (diffField)
		{
			for (auto K : diffField->Tensors())
				CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreReconstructStiffnessMatrix(*K, reconstructionBasis);
		}
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
#ifdef ENABLE_3D
			// - Tetrahedral element
			Tetrahedron::InitReferenceShape()->ComputeAndStoreMassMatrix((FunctionalBasis<3>*)cellBasis);
			Tetrahedron::InitReferenceShape()->ComputeAndStoreMassMatrix((FunctionalBasis<3>*)reconstructionBasis);
			Tetrahedron::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix((FunctionalBasis<3>*)cellBasis, (FunctionalBasis<3>*)reconstructionBasis);
			// - Triangular face
			Triangle::InitReferenceShape()->ComputeAndStoreMassMatrix((FunctionalBasis<2>*)faceBasis);
#endif // ENABLE_3D
		}
	}

	//------------------------------//
	//           Init HHO           //
	//------------------------------//
	void InitHHO(bool assembleLocalMatrices = true)
	{
		InitHHO_Faces();
		InitHHO_Elements(assembleLocalMatrices);

		SkeletonSpace    = HHOSkeletonSpace   (_mesh, HHO, _hhoFaces);
		ReconstructSpace = HHOReconstructSpace(_mesh, HHO, _hhoElements);
		BoundarySpace    = HHOBoundarySpace   (_mesh, HHO, _hhoFaces);
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
	void InitHHO_Elements(bool assembleLocalMatrices = true)
	{
		if (!_hhoElements.empty())
			Utils::Warning("The HHO elements of this problem have already been initialized.");

		_hhoElements = vector<Diff_HHOElement<Dim>>(this->_mesh->Elements.size());

		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, assembleLocalMatrices](Element<Dim>* e)
			{
				_hhoElements[e->Number].MeshElement = e;

				for (Face<Dim>* f : e->Faces)
					_hhoElements[e->Number].Faces.push_back(&this->_hhoFaces[f->Number]);

				_hhoElements[e->Number].InitHHO(HHO, assembleLocalMatrices);
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

	void DeleteHHOElements()
	{
		this->_hhoElements.clear();
	}
	void DeleteHHOFaces()
	{
		this->_hhoFaces.clear();
	}

	//------------------------------//
	//   Assemble right-hand side   //
	//------------------------------//

	Vector AssembleSourceTerm(DomFunction sourceFunction)
	{
		Vector b_source = Vector(HHO->nTotalCellUnknowns);
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &b_source, &sourceFunction](Element<Dim>* e)
			{
				Diff_HHOElement<Dim>* element = HHOElement(e);
				for (BasisFunction<Dim>* cellPhi : element->CellBasis->LocalFunctions)
				{
					BigNumber i = DOFNumber(element, cellPhi);
					b_source(i) = element->SourceTerm(cellPhi, sourceFunction);
				}
			}
		);
		return b_source;
	}

	Vector AssembleSourceTerm(const Vector& sourceFuncCoeffs)
	{
		ElementParallelLoop<Dim> parallelLoop(this->_mesh->Elements);
		Vector b_source = Vector(HHO->nTotalCellUnknowns);

		if (sourceFuncCoeffs.rows() == HHO->nTotalCellUnknowns) // degree k
		{
			parallelLoop.Execute([this, &sourceFuncCoeffs, &b_source](Element<Dim>* e)
				{
					Diff_HHOElement<Dim>* element = this->HHOElement(e);
					b_source.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns) = element->ApplyCellMassMatrix(sourceFuncCoeffs.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns));
				}
			);
		}
		else if (sourceFuncCoeffs.rows() == HHO->nTotalReconstructUnknowns) // degree k+1
		{
			parallelLoop.Execute([this, &sourceFuncCoeffs, &b_source](Element<Dim>* e)
				{
					Diff_HHOElement<Dim>* element = this->HHOElement(e);
					b_source.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns) = element->ApplyCellReconstructMassMatrix(sourceFuncCoeffs.segment(e->Number * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns));
				}
			);
		}
		else
			Utils::FatalError("the argument sourceFuncCoeffs does not have a correct size");

		return b_source;
	}

	// Solution on the Dirichlet faces
	Vector AssembleDirichletTerm(DomFunction dirichletFunction)
	{
		Vector x_dF = Vector(HHO->nDirichletCoeffs);
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->DirichletFaces, [this, &x_dF, &dirichletFunction](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = FirstDOFGlobalNumber(face) - HHO->nTotalHybridUnknowns;
				x_dF.segment(i, HHO->nFaceUnknowns) = face->ProjectOnBasis(dirichletFunction);
			}
		);
		return x_dF;
	}

	/*Vector AssembleDirichletTerm(const Vector& x_dF)
	{
		if (x_dF.rows() == HHO->nDirichletFaces * HHO->nFaceUnknowns)
			return x_dF;
		else if ()
	}*/

	Vector AssembleDirichletTermFromReconstructedBoundaryElem(const Vector& reconstructedElem)
	{
		assert(reconstructedElem.rows() == _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);

		Vector x_dF = Vector(HHO->nDirichletCoeffs);
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->DirichletFaces, [this, &x_dF, &reconstructedElem](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				Diff_HHOElement<Dim>* elem = HHOElement(f->Element1);
				BigNumber i = FirstDOFGlobalNumber(face) - HHO->nTotalHybridUnknowns;
				int boundaryElemNumber = _mesh->BoundaryElementNumber(elem->MeshElement);
				x_dF.segment(i, HHO->nFaceUnknowns) = face->Trace(elem->MeshElement, elem->ReconstructionBasis)* reconstructedElem.segment(boundaryElemNumber * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns);
			}
		);
		return x_dF;
	}
	
	Vector ReconstructAndAssembleDirichletTerm(const Vector& boundaryCoeffs)
	{
		assert(boundaryCoeffs.rows() == HHO->nDirichletCoeffs);

		Vector x_dF = Vector(HHO->nDirichletCoeffs);
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->DirichletFaces, [this, &x_dF, &boundaryCoeffs](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = FirstDOFGlobalNumber(face) - HHO->nTotalHybridUnknowns;
				Element<Dim>* e = f->Element1;
				Diff_HHOElement<Dim>* elem = HHOElement(e);
				// Reconstruct k+1
				Vector faceCoeffs = Vector::Zero(e->Faces.size() * HHO->nFaceUnknowns);
				for (Face<Dim>* ef : e->Faces)
				{
					if (ef->IsDomainBoundary)
						faceCoeffs.segment(e->LocalNumberOf(ef) * HHO->nFaceUnknowns, HHO->nFaceUnknowns) = boundaryCoeffs.segment((ef->Number - HHO->nInteriorFaces) * HHO->nFaceUnknowns, HHO->nFaceUnknowns);
				}
				Vector reconstruct = elem->ReconstructFromFaces(faceCoeffs);
				// Trace of degree k
				x_dF.segment(i, HHO->nFaceUnknowns) = face->Trace(e, elem->ReconstructionBasis) * reconstruct;
			}
		);
		return x_dF;
	}

	//         [    0    ]   <---- interior faces
	// returns [---------]
	//         [ b_neumF ]   <---- Neumann faces
	// 
	// where b_neumF = [ (neumannFunc|phi_i)_F ]
	Vector AssembleNeumannTerm(DomFunction neumannFunction)
	{
		Vector b_neumann = Vector::Zero(HHO->nTotalFaceUnknowns);
		//this->B_neumF = Vector(HHO->nNeumannUnknowns);
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->NeumannFaces, [this, &b_neumann, &neumannFunction](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = FirstDOFGlobalNumber(face) - HHO->nTotalCellUnknowns;
				b_neumann.segment(i, HHO->nFaceUnknowns) = face->InnerProductWithBasis(neumannFunction);
				//BigNumber i = FirstDOFGlobalNumber(face) - HHO->nTotalCellUnknowns - HHO->nInteriorFaces*HHO->nFaceUnknowns;
				//this->B_neumF.segment(i, HHO->nFaceUnknowns) = face->InnerProductWithBasis(neumannFunction);
			}
		);
		return b_neumann;
	}

	//         [    0    ]   <---- interior faces
	// returns [---------]
	//         [ b_neumF ]   <---- Neumann faces
	// 
	// where b_neumF = [ M * coeffs ]
	Vector AssembleNeumannTerm(const Vector& neumannFuncCoeffs)
	{
		assert(neumannFuncCoeffs.rows() == HHO->nNeumannUnknowns);
		assert(neumannFuncCoeffs.rows() == this->_mesh->NeumannFaces.size() * HHO->nFaceUnknowns);

		Vector b_neumann = Vector::Zero(HHO->nTotalFaceUnknowns);
		ParallelLoop<Face<Dim>*>::Execute(this->_mesh->NeumannFaces, [this, &b_neumann, &neumannFuncCoeffs](Face<Dim>* f)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				BigNumber i = face->Number() * HHO->nFaceUnknowns;
				assert(i >= HHO->nInteriorFaces * HHO->nFaceUnknowns);
				BigNumber j = (face->Number() - HHO->nInteriorFaces) * HHO->nFaceUnknowns;
				b_neumann.segment(i, HHO->nFaceUnknowns) = face->ApplyMassMatrix(neumannFuncCoeffs.segment(j, HHO->nFaceUnknowns));
			}
		);
		return b_neumann;
	}

	//-------------------------------------------------------------------------------------------------//
	// After solving the faces, construction of the higher-order approximation using the reconstructor //
	//-------------------------------------------------------------------------------------------------//

	Vector HybridCoeffsBySolvingCellUnknowns(const Vector& faceUnknowns)
	{
		return HybridCoeffsBySolvingCellUnknowns(faceUnknowns, this->x_dF, this->B_T);
	}

	Vector HybridCoeffsBySolvingCellUnknowns(const Vector& faceUnknowns, const Vector& dirichletCoeffs, const Vector& b_T)
	{
		assert(faceUnknowns.rows() == HHO->nTotalFaceUnknowns);
		assert(dirichletCoeffs.rows() == HHO->nDirichletCoeffs);

		Vector hybridCoeffs = Vector(HHO->nTotalHybridCoeffs);
		hybridCoeffs.head(HHO->nTotalCellUnknowns) = SolveCellUnknowns(faceUnknowns, b_T);
		hybridCoeffs.segment(HHO->nTotalCellUnknowns, HHO->nTotalFaceUnknowns) = faceUnknowns;
		hybridCoeffs.tail(HHO->nDirichletCoeffs) = dirichletCoeffs; // Dirichlet boundary conditions

		return hybridCoeffs;
	}

	Vector HybridCoeffsByAddingDirichletCoeffs(const Vector& hybridUnknowns)
	{
		assert(hybridUnknowns.rows() == HHO->nTotalHybridUnknowns);

		Vector hybridCoeffs = Vector(HHO->nTotalHybridCoeffs);
		hybridCoeffs.head(HHO->nTotalHybridUnknowns) = hybridUnknowns;
		hybridCoeffs.tail(HHO->nDirichletCoeffs) = this->x_dF; // Dirichlet boundary conditions

		return hybridCoeffs;
	}

	Vector ReconstructHigherOrderApproximationFromFaceCoeffs(const Vector& faceUnknowns)
	{
		Vector hybridCoeffs = HybridCoeffsBySolvingCellUnknowns(faceUnknowns);
		return ReconstructHigherOrderApproximationFromHybridCoeffs(hybridCoeffs);
	}

	Vector ReconstructHigherOrderApproximationFromFaceCoeffs(const Vector& faceUnknowns, const Vector& dirichletCoeffs, const Vector& b_T)
	{
		Vector hybridCoeffs = HybridCoeffsBySolvingCellUnknowns(faceUnknowns, dirichletCoeffs, b_T);
		return ReconstructHigherOrderApproximationFromHybridCoeffs(hybridCoeffs);
	}

	Vector ReconstructHigherOrderApproximationFromHybridCoeffs(const Vector& hybridCoeffs)
	{
		assert(hybridCoeffs.rows() == HHO->nTotalHybridCoeffs);
		Vector reconstruction(HHO->nElements * HHO->nReconstructUnknowns);
		ParallelLoop<Element<Dim>*>::Execute(this->_mesh->Elements, [this, &reconstruction, &hybridCoeffs](Element<Dim>* e)
			{
				Diff_HHOElement<Dim>* elem = HHOElement(e);

				Vector localHybrid(HHO->nCellUnknowns + HHO->nFaceUnknowns * elem->Faces.size());
				localHybrid.head(HHO->nCellUnknowns) = hybridCoeffs.segment(FirstDOFGlobalNumber(elem), HHO->nCellUnknowns);
				for (auto face : elem->Faces)
					localHybrid.segment(elem->FirstDOFNumber(face), HHO->nFaceUnknowns) = hybridCoeffs.segment(FirstDOFGlobalNumber(face), HHO->nFaceUnknowns);

				Vector localReconstruction = elem->Reconstruct(localHybrid);
				reconstruction.segment(elem->Number() * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns) = localReconstruction;
			});
		return reconstruction;
	}

	Vector ReconstructHigherOrderOnBoundaryOnly(const Vector& faceUnknowns, const Vector& dirichletCoeffs, const Vector& b_T)
	{
		Vector v = b_T - A_T_ndF * faceUnknowns;

		Vector faceCoeffs(HHO->nTotalFaceUnknowns + HHO->nDirichletCoeffs);
		faceCoeffs.head(HHO->nTotalFaceUnknowns) = faceUnknowns;
		faceCoeffs.tail(HHO->nDirichletCoeffs)   = dirichletCoeffs;

		Vector reconstruction(_mesh->NBoundaryElements() * HHO->nReconstructUnknowns);

		ElementParallelLoop<Dim> parallelLoop(this->_mesh->Elements);
		parallelLoop.Execute([this, &faceCoeffs, &v, &reconstruction](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				if (e->IsOnBoundary())
				{
					Diff_HHOElement<Dim>* elem = HHOElement(e);

					Vector localHybrid(HHO->nCellUnknowns + HHO->nFaceUnknowns * e->Faces.size());
					localHybrid.head(HHO->nCellUnknowns) = elem->AttSolver.solve(v.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns));
					for (auto face : elem->Faces)
						localHybrid.segment(elem->FirstDOFNumber(face), HHO->nFaceUnknowns) = faceCoeffs.segment(face->Number() * HHO->nFaceUnknowns, HHO->nFaceUnknowns);

					int boundaryElemNumber = _mesh->BoundaryElementNumber(e);
					reconstruction.segment(boundaryElemNumber * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns) = elem->Reconstruct(localHybrid);
				}
			});
		return reconstruction;
	}

	Vector SolveCellUnknownsOnBoundaryOnly(const Vector& faceUnknowns, const Vector& b_T)
	{
		Vector v = b_T - A_T_ndF * faceUnknowns;

		Vector cellUnknowns(_mesh->NBoundaryElements() * HHO->nCellUnknowns);

		ElementParallelLoop<Dim> parallelLoop(this->_mesh->Elements);
		parallelLoop.Execute([this, &v, &cellUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				if (e->IsOnBoundary())
				{
					Diff_HHOElement<Dim>* elem = HHOElement(e);

					int boundaryElemNumber = _mesh->BoundaryElementNumber(e);
					cellUnknowns.segment(boundaryElemNumber * HHO->nCellUnknowns, HHO->nCellUnknowns) = elem->AttSolver.solve(v.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns));
				}
			});
		return cellUnknowns;
	}

	Vector SolveCellUnknowns(const Vector& faceUnknowns, const Vector& b_T)
	{
		assert(faceUnknowns.rows() == HHO->nTotalFaceUnknowns);
		assert(b_T.rows() == HHO->nTotalCellUnknowns);

		return Solve_A_T_T(b_T - A_T_ndF * faceUnknowns);
	}

	// Deprecated!
	// Reassembles B_T w.r.t. to a discrete source function
	void ChangeSourceFunction(const Vector& sourceFuncCoeffs)
	{
		assert(this->x_dF.rows() == 0 || this->x_dF.isZero(0));

		this->B_T = std::move(AssembleSourceTerm(sourceFuncCoeffs));
		
		// Elimination of Dirichlet unknowns (here 0 (so far), but non-homogeneous Dirichlet BC must be managed!)
		/*
		SparseMatrix A_ndF_dF = A_F_F.topRightCorner(HHO->nTotalFaceUnknowns, HHO->nDirichletCoeffs);
		SparseMatrix A_T_dF   = A_T_F.topRightCorner(HHO->nTotalCellUnknowns, HHO->nDirichletCoeffs);

		this->B_T   -= A_T_dF   * this->x_dF;
		this->B_ndF -= A_ndF_dF * this->x_dF;
		*/
	}

	//-------------------------//
	//     Right-hand side     //
	//-------------------------//

	Vector ComputeB_T(const Vector& b_source, const Vector& x_dF)
	{
		assert(b_source.rows() == HHO->nTotalCellUnknowns);
		assert(x_dF.rows() == HHO->nDirichletCoeffs);

		return b_source - A_T_dF * x_dF;
	}

	Vector ComputeB_T_zeroSource(const Vector& x_dF)
	{
		assert(x_dF.rows() == HHO->nDirichletCoeffs);

		return -A_T_dF * x_dF;
	}

	Vector ComputeB_ndF(const Vector& b_neumann, const Vector& x_dF)
	{
		assert(x_dF.rows() == HHO->nDirichletCoeffs);

		return b_neumann - A_ndF_dF * x_dF;
	}

	Vector ComputeB_ndF_noNeumann(const Vector& x_dF)
	{
		assert(x_dF.rows() == HHO->nDirichletCoeffs);

		return -A_ndF_dF * x_dF;
	}

	// Deprecated
	Vector& SetCondensedRHS()
	{
		this->b = this->B_ndF -this->A_T_ndF.transpose() * Solve_A_T_T(this->B_T);
		return this->b;
	}

	Vector CondensedRHS(const Vector& b_T, const Vector& b_ndF)
	{
		assert(b_T.rows() == HHO->nTotalCellUnknowns);
		assert(b_ndF.rows() == HHO->nTotalFaceUnknowns);

		return b_ndF - this->A_T_ndF.transpose() * Solve_A_T_T(b_T);
	}

	Vector CondensedRHS_noNeumannZeroDirichlet(const Vector& b_T)
	{
		assert(b_T.rows() == HHO->nTotalCellUnknowns);

		return -this->A_T_ndF.transpose() * Solve_A_T_T(b_T);
	}

	Vector CondensedRHS_noDirichletZeroNeumann(const Vector& b_T)
	{
		assert(b_T.rows() == HHO->nTotalCellUnknowns);

		return -this->A_T_ndF.transpose() * Solve_A_T_T(b_T);
	}

	// Returns the matrix that extract the normal derivative on the boundary from a reconstructed vector restricted to the boundary elements.
	// The result is expressed w.r.t. the face unknowns. 
	// It is possible because the derivative a reconstructed polynomial is (k+1)-1, i.e. the degree of the face unknowns.
	SparseMatrix NormalDerivativeMatrix()
	{
		FaceParallelLoop<Dim> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nFaceUnknowns);

		parallelLoop.Execute([this](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				int i = f->Number - HHO->nInteriorFaces;
				Diff_HHOFace<Dim>* face = HHOFace(f);

				int j = _mesh->BoundaryElementNumber(f->Element1);
				Diff_HHOElement<Dim>* elem = this->HHOElement(f->Element1);

				chunk->Results.Coeffs.Add(i * HHO->nFaceUnknowns, j * HHO->nReconstructUnknowns, face->NormalDerivative(elem->MeshElement, elem->ReconstructionBasis));
			});

		SparseMatrix normalDerivative(HHO->nBoundaryFaces * HHO->nFaceUnknowns, _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		parallelLoop.Fill(normalDerivative);
		return normalDerivative;
	}

	// Computes P^T*A*P*x where A is the stiffness matrix and the argument is P*x
	SparseMatrix PTranspose_Stiff_boundary()
	{
		FaceParallelLoop<Dim> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nFaceUnknowns);

		parallelLoop.Execute([this](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				int i = f->Number - HHO->nInteriorFaces;
				Diff_HHOFace<Dim>* face = HHOFace(f);

				int j = _mesh->BoundaryElementNumber(f->Element1);
				Diff_HHOElement<Dim>* elem = this->HHOElement(f->Element1);

				//DenseMatrix P = elem->P.middleCols(HHO->nCellUnknowns + elem->LocalNumberOf(face) * HHO->nFaceUnknowns, HHO->nFaceUnknowns);
				DenseMatrix P = elem->ReconstructionFromFacesMatrix().middleCols(elem->LocalNumberOf(face) * HHO->nFaceUnknowns, HHO->nFaceUnknowns);
				DenseMatrix m = P.transpose() * elem->MeshElement->IntegralGradGradMatrix(elem->ReconstructionBasis);

				chunk->Results.Coeffs.Add(i * HHO->nFaceUnknowns, j * HHO->nReconstructUnknowns, m);
			});

		SparseMatrix mat(HHO->nBoundaryFaces * HHO->nFaceUnknowns, _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		parallelLoop.Fill(mat);
		return mat;
	}

	// Computes P^T*M*P*x where the argument is P*x
	SparseMatrix PTranspose_Mass_boundary()
	{
		FaceParallelLoop<Dim> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nFaceUnknowns);

		parallelLoop.Execute([this](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				int i = f->Number - HHO->nInteriorFaces;
				Diff_HHOFace<Dim>* face = HHOFace(f);

				int j = f->Element1->Number;
				Diff_HHOElement<Dim>* elem = this->HHOElement(f->Element1);

				//DenseMatrix P = elem->P.middleCols(HHO->nCellUnknowns + elem->LocalNumberOf(face) * HHO->nFaceUnknowns, HHO->nFaceUnknowns);
				DenseMatrix P = elem->ReconstructionFromFacesMatrix().middleCols(elem->LocalNumberOf(face) * HHO->nFaceUnknowns, HHO->nFaceUnknowns);
				DenseMatrix m = P.transpose() * elem->MassMatrix(elem->ReconstructionBasis);

				chunk->Results.Coeffs.Add(i * HHO->nFaceUnknowns, j * HHO->nReconstructUnknowns, m);
			});

		SparseMatrix mat(HHO->nBoundaryFaces * HHO->nFaceUnknowns, HHO->nTotalReconstructUnknowns);
		parallelLoop.Fill(mat);
		return mat;
	}

	SparseMatrix PTranspose()
	{
		FaceParallelLoop<Dim> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nFaceUnknowns);

		parallelLoop.Execute([this](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				int i = f->Number - HHO->nInteriorFaces;
				Diff_HHOFace<Dim>* face = HHOFace(f);

				//int j = _mesh->BoundaryElementNumber(f->Element1);
				Diff_HHOElement<Dim>* elem = this->HHOElement(f->Element1);

				DenseMatrix P = elem->ReconstructionFromFacesMatrix().middleCols(elem->LocalNumberOf(face) * HHO->nFaceUnknowns, HHO->nFaceUnknowns);

				chunk->Results.Coeffs.Add(i * HHO->nFaceUnknowns, elem->Number() * HHO->nReconstructUnknowns, P.transpose());
			});

		SparseMatrix mat(HHO->nBoundaryFaces * HHO->nFaceUnknowns, HHO->nTotalReconstructUnknowns);
		parallelLoop.Fill(mat);
		return mat;
	}

	SparseMatrix PTransposeOnBoundary()
	{
		FaceParallelLoop<Dim> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nFaceUnknowns);

		parallelLoop.Execute([this](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				int i = f->Number - HHO->nInteriorFaces;
				Diff_HHOFace<Dim>* face = HHOFace(f);

				int j = _mesh->BoundaryElementNumber(f->Element1);
				Diff_HHOElement<Dim>* elem = this->HHOElement(f->Element1);

				DenseMatrix P = elem->ReconstructionFromFacesMatrix().middleCols(elem->LocalNumberOf(face) * HHO->nFaceUnknowns, HHO->nFaceUnknowns);

				chunk->Results.Coeffs.Add(i * HHO->nFaceUnknowns, j * HHO->nReconstructUnknowns, P.transpose());
			});

		SparseMatrix mat(HHO->nBoundaryFaces * HHO->nFaceUnknowns, _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		parallelLoop.Fill(mat);
		return mat;
	}

	SparseMatrix SolveCellUknTransposeOnBoundary()
	{
		FaceParallelLoop<Dim> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nCellUnknowns * HHO->nFaceUnknowns);

		parallelLoop.Execute([this](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				int i = f->Number - HHO->nInteriorFaces;
				Diff_HHOFace<Dim>* face = HHOFace(f);

				int j = _mesh->BoundaryElementNumber(f->Element1);
				Diff_HHOElement<Dim>* elem = this->HHOElement(f->Element1);

				DenseMatrix S = elem->SolveCellUnknownsMatrix().middleCols(elem->LocalNumberOf(face) * HHO->nFaceUnknowns, HHO->nFaceUnknowns);

				chunk->Results.Coeffs.Add(i * HHO->nFaceUnknowns, j * HHO->nCellUnknowns, S.transpose());
			});

		SparseMatrix mat(HHO->nBoundaryFaces * HHO->nFaceUnknowns, _mesh->NBoundaryElements() * HHO->nCellUnknowns);
		parallelLoop.Fill(mat);
		return mat;
	}

	SparseMatrix LowerOrderTraceMatrixOnBoundary()
	{
		FaceParallelLoop<Dim> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nFaceUnknowns);

		parallelLoop.Execute([this](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				int i = f->Number - HHO->nInteriorFaces;
				Diff_HHOFace<Dim>* face = HHOFace(f);

				int j = _mesh->BoundaryElementNumber(f->Element1);
				Diff_HHOElement<Dim>* elem = this->HHOElement(f->Element1);

				chunk->Results.Coeffs.Add(i * HHO->nFaceUnknowns, j * HHO->nReconstructUnknowns, face->Trace(f->Element1, elem->ReconstructionBasis));
			});

		SparseMatrix mat(HHO->nBoundaryFaces * HHO->nFaceUnknowns, _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		parallelLoop.Fill(mat);
		return mat;
	}

	// Stiffness matrix (grad, grad) of degree k+1
	SparseMatrix ReconstructStiffnessMatrix()
	{
		ElementParallelLoop<Dim> parallelLoop(_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nReconstructUnknowns);

		parallelLoop.Execute([this](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* elem = this->HHOElement(e);

				chunk->Results.Coeffs.Add(e->Number * HHO->nReconstructUnknowns, e->Number * HHO->nReconstructUnknowns, e->IntegralGradGradMatrix(elem->ReconstructionBasis));
			});

		SparseMatrix stiff(HHO->nTotalReconstructUnknowns, HHO->nTotalReconstructUnknowns);
		parallelLoop.Fill(stiff);
		return stiff;
	}

	// Stiffness matrix (grad, grad) of degree k+1, on boundary elements only
	SparseMatrix ReconstructStiffnessMatrixOnBoundaryElements()
	{
		ElementParallelLoop<Dim> parallelLoop(_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nReconstructUnknowns);

		parallelLoop.Execute([this](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				if (e->IsOnBoundary())
				{
					int i = _mesh->BoundaryElementNumber(e);
					Diff_HHOElement<Dim>* elem = this->HHOElement(e);

					chunk->Results.Coeffs.Add(i * HHO->nReconstructUnknowns, i * HHO->nReconstructUnknowns, elem->MeshElement->IntegralGradGradMatrix(elem->ReconstructionBasis));
				}
			});

		SparseMatrix stiff(_mesh->NBoundaryElements() * HHO->nReconstructUnknowns, _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		parallelLoop.Fill(stiff);

		//cout << DenseMatrix(stiff) << endl;

		return stiff;
	}

	// Stiffness matrix (grad, grad) of degree k, on boundary elements only
	SparseMatrix CellStiffnessMatrixOnBoundaryElements()
	{
		ElementParallelLoop<Dim> parallelLoop(_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nCellUnknowns * HHO->nCellUnknowns);

		parallelLoop.Execute([this](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				if (e->IsOnBoundary())
				{
					int i = _mesh->BoundaryElementNumber(e);
					Diff_HHOElement<Dim>* elem = this->HHOElement(e);

					chunk->Results.Coeffs.Add(i * HHO->nCellUnknowns, i * HHO->nCellUnknowns, elem->MeshElement->IntegralGradGradMatrix(elem->CellBasis));
				}
			});

		SparseMatrix stiff(_mesh->NBoundaryElements() * HHO->nCellUnknowns, _mesh->NBoundaryElements() * HHO->nCellUnknowns);
		parallelLoop.Fill(stiff);
		return stiff;
	}

	// Mass matrix of degree k+1
	SparseMatrix ReconstructMassMatrix()
	{
		SparseMatrix mass(HHO->nTotalReconstructUnknowns, HHO->nTotalReconstructUnknowns);
		if (HHO->OrthonormalizeElemBases())
			mass.setIdentity();
		else
		{
			ElementParallelLoop<Dim> parallelLoop(_mesh->Elements);
			parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nReconstructUnknowns);

			parallelLoop.Execute([this](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
				{
					Diff_HHOElement<Dim>* elem = this->HHOElement(e);
					chunk->Results.Coeffs.Add(e->Number * HHO->nReconstructUnknowns, e->Number * HHO->nReconstructUnknowns, elem->MassMatrix(elem->ReconstructionBasis));
				});

			parallelLoop.Fill(mass);
		}
		return mass;
	}

	// Mass matrix of degree k+1, on boundary elements only
	SparseMatrix ReconstructMassMatrixOnBoundaryElements()
	{
		SparseMatrix mass(_mesh->NBoundaryElements() * HHO->nReconstructUnknowns, _mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
		if (HHO->OrthonormalizeElemBases())
			mass.setIdentity();
		else
		{
			ElementParallelLoop<Dim> parallelLoop(_mesh->Elements);
			parallelLoop.ReserveChunkCoeffsSize(HHO->nReconstructUnknowns * HHO->nReconstructUnknowns);

			parallelLoop.Execute([this](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
				{
					if (e->IsOnBoundary())
					{
						int i = _mesh->BoundaryElementNumber(e);
						Diff_HHOElement<Dim>* elem = this->HHOElement(e);

						chunk->Results.Coeffs.Add(i * HHO->nReconstructUnknowns, i * HHO->nReconstructUnknowns, elem->MassMatrix(elem->ReconstructionBasis));
					}
				});
			parallelLoop.Fill(mass);
		}
		//cout << DenseMatrix(mass) << endl;
		return mass;
	}

	// Mass matrix of degree k, on boundary elements only
	SparseMatrix CellMassMatrixOnBoundaryElements()
	{
		SparseMatrix mass(_mesh->NBoundaryElements() * HHO->nCellUnknowns, _mesh->NBoundaryElements() * HHO->nCellUnknowns);
		if (HHO->OrthonormalizeElemBases())
			mass.setIdentity();
		else
		{
			ElementParallelLoop<Dim> parallelLoop(_mesh->Elements);
			parallelLoop.ReserveChunkCoeffsSize(HHO->nCellUnknowns * HHO->nCellUnknowns);

			parallelLoop.Execute([this](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
				{
					if (e->IsOnBoundary())
					{
						int i = _mesh->BoundaryElementNumber(e);
						Diff_HHOElement<Dim>* elem = this->HHOElement(e);

						chunk->Results.Coeffs.Add(i * HHO->nCellUnknowns, i * HHO->nCellUnknowns, elem->MassMatrix(elem->CellBasis));
					}
				});

			parallelLoop.Fill(mass);
		}
		return mass;
	}

	Vector ExtractElemBoundary(const Vector v)
	{
		if (v.rows() == HHO->nTotalReconstructUnknowns)
		{
			Vector boundary(_mesh->NBoundaryElements() * HHO->nReconstructUnknowns);
			for (Element<Dim>* e : _mesh->Elements)
			{
				if (e->IsOnBoundary())
				{
					int i = _mesh->BoundaryElementNumber(e);
					boundary.segment(i * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns) = v.segment(e->Number * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns);
				}
			}
			return boundary;
		}
		else if (v.rows() == HHO->nTotalCellUnknowns)
		{
			Vector boundary(_mesh->NBoundaryElements() * HHO->nCellUnknowns);
			for (Element<Dim>* e : _mesh->Elements)
			{
				if (e->IsOnBoundary())
				{
					int i = _mesh->BoundaryElementNumber(e);
					boundary.segment(i * HHO->nCellUnknowns, HHO->nCellUnknowns) = v.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns);
				}
			}
			return boundary;
		}
		Utils::FatalError("Not implemented");
		return Vector();
	}

	Vector SolveFaceMassMatrixOnBoundary(const Vector& v)
	{
		Vector res(v.rows());
		FaceParallelLoop<Dim> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.Execute([this, &v, &res](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOFace<Dim>* face = HHOFace(f);
				int i = f->Number - HHO->nInteriorFaces;

				res.segment(i * HHO->nFaceUnknowns, HHO->nFaceUnknowns) = face->SolveMassMatrix(v.segment(i * HHO->nFaceUnknowns, HHO->nFaceUnknowns));
			});
		return res;
	}

	//---------------------------------------//
	//                Exports                //
	//---------------------------------------//

public:
	void ExportReconstructedVectorToGMSH(const Vector& reconstructedVector, const ExportModule& out, string name, double tolerance = 1e-3, int maxRefinements = 6)
	{
		if (HHO->OrthogonalizeElemBases())
			Utils::Error("The export to GMSH has not been implemented when the bases are orthonormalized against each element.");
		this->_mesh->ExportToGMSH_Elements(this->HHO->ReconstructionBasis, reconstructedVector, out.GetFilePathPrefix(), name, tolerance, maxRefinements);
	}

	void ExportSolutionToGMSH(const Vector& reconstructedSolution, const ExportModule& out, double tolerance = 1e-3, int maxRefinements = 6)
	{
		ExportReconstructedVectorToGMSH(reconstructedSolution, out, "potential", tolerance, maxRefinements);
	}

	void ExportErrorToGMSH(const Vector& faceCoeffs, const ExportModule& out)
	{
		if (HHO->OrthogonalizeElemBases())
			Utils::Error("The export to GMSH has not been implemented when the bases are orthonormalized against each element.");
		Vector cellCoeffs = Solve_A_T_T(B_T - A_T_ndF * faceCoeffs);
		this->_mesh->ExportToGMSH_Elements(this->HHO->CellBasis, cellCoeffs, out.GetFilePathPrefix(), "error");
	}

	//---------------------------------------//
	//              Solve_A_T_T              //
	//---------------------------------------//

	SparseMatrix Solve_A_T_T(const SparseMatrix& A_T_ndF)
	{
		ElementParallelLoop<Dim> parallelLoop(this->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nCellUnknowns * HHO->nCellUnknowns);
		parallelLoop.Execute([this, &A_T_ndF](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* element = HHOElement(e);
				for (auto f : element->MeshElement->Faces)
				{
					if (f->HasDirichletBC())
						continue;
					DenseMatrix block_A_T_ndF = A_T_ndF.block(e->Number * HHO->nCellUnknowns, f->Number * HHO->nFaceUnknowns, HHO->nCellUnknowns, HHO->nFaceUnknowns);
					chunk->Results.Coeffs.Add(e->Number * HHO->nCellUnknowns, f->Number * HHO->nFaceUnknowns, element->AttSolver.solve(block_A_T_ndF));
				}
			});
		SparseMatrix result = SparseMatrix(HHO->nTotalCellUnknowns, HHO->nTotalFaceUnknowns);
		parallelLoop.Fill(result);
		return result;
	}

	Vector Solve_A_T_T(const Vector& b_T)
	{
		Vector result(b_T.rows());
		ElementParallelLoop<Dim> parallelLoop(this->_mesh->Elements);
		parallelLoop.Execute([this, &b_T, &result](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* element = HHOElement(e);
				result.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns) = element->AttSolver.solve(b_T.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns));
			});
		return result;
	}

	//--------------------------------------------//
	//              Global integrals              //
	//--------------------------------------------//

	/*double IntegralFromFaceDoFs(Vector dofs, int polyDegree)
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
	}*/

	double IntegralOverDomain(DomFunction func)
	{
		struct ChunkResult { double total = 0; };

		ParallelLoop<Element<Dim>*, ChunkResult> parallelLoop(_mesh->Elements);
		parallelLoop.Execute([this, func](Element<Dim>* e, ParallelChunk<ChunkResult>* chunk)
			{
				chunk->Results.total += e->Integral(func);
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
	}

	double IntegralOverBoundary(DomFunction func)
	{
		struct ChunkResult { double total = 0; };

		ParallelLoop<Face<Dim>*, ChunkResult> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.Execute([this, func](Face<Dim>* f, ParallelChunk<ChunkResult>* chunk)
			{
				chunk->Results.total += f->Integral(func);
			});

		double total = 0;
		parallelLoop.AggregateChunkResults([&total](ChunkResult chunkResult)
			{
				total += chunkResult.total;
			});
		return total;
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

