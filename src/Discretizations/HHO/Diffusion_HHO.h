#pragma once
#include "../../Mesh/Mesh.h"
#include "../../Utils/Utils.h"
#ifdef ENABLE_3D
#include "../../Geometry/3D/Tetrahedron.h"
#endif // ENABLE_3D
#include "../../Geometry/CartesianShape.h"
#include "../../Geometry/2D/Triangle.h"
#include "../../Geometry/2D/Quadrilateral.h"
#include "DiscreteSpaces/HHOCellSpace.h"
#include "DiscreteSpaces/HHOSkeletonSpace.h"
#include "DiscreteSpaces/HHOReconstructSpace.h"
#include "DiscreteSpaces/HHOBoundarySpace.h"
#include "DiscreteSpaces/HHODirichletSpace.h"
#include "DiscreteSpaces/HHONeumannSpace.h"
#include "DiscreteSpaces/HHONonDirichletFaceSpace.h"
#include "Diffusion_HHOMatrix.h"
#include "../../TestCases/Diffusion/DiffusionTestCase.h"
#include "../../Utils/ExportModule.h"
using namespace std;

template <int Dim>
class Diffusion_HHO
{
	template <int Dim2>
	friend class NeighbourhoodDiffusion_HHO;

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

	HHOCellSpace<Dim>             CellSpace;
	HHOReconstructSpace<Dim>      ReconstructSpace;
	HHOSkeletonSpace<Dim>		  SkeletonSpace;
	HHOBoundarySpace<Dim>         BoundarySpace;
	HHODirichletSpace<Dim>        DirichletSpace;
	HHONeumannSpace<Dim>          NeumannSpace;
	HHONonDirichletFaceSpace<Dim> NonDirichletFaceSpace;

	// Matrix parts:
	// Used in reconstruction: A_T_ndF
	// Used in right-hand side: A_T_dF, A_ndF_dF
	SparseMatrix A_T_T, A_T_ndF,   A_T_dF;
	SparseMatrix        A_ndF_ndF, A_ndF_dF;
	SparseMatrix                   A_dF_dF;
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
		FunctionalBasis<Dim>* reconstructionBasis = HHO->ReconstructionBasis->CreateLowerDegreeBasis(newReconstructDegree);
		FunctionalBasis<Dim>* cellBasis = HHO->CellBasis->CreateLowerDegreeBasis(newCellDegree);
		FunctionalBasis<Dim-1>* faceBasis = HHO->FaceBasis->CreateLowerDegreeBasis(faceDegree);

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
		FunctionalBasis<Dim>* reconstructionBasis = HHO->ReconstructionBasis->CreateLowerDegreeBasis(newReconstructDegree);
		FunctionalBasis<Dim>* cellBasis = HHO->CellBasis->CreateLowerDegreeBasis(newCellDegree);
		FunctionalBasis<Dim - 1>* faceBasis = HHO->FaceBasis->CreateLowerDegreeBasis(faceDegree);

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

		parallelLoop.InitChunks([this, actions, nCoeffs_A_T, nCoeffs_A_T_T, nCoeffs_A_T_F, nCoeffs_A_F_F](ParallelChunk<AssemblyResult>* chunk)
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

		if (this->TestCase->BC.Type != PbBoundaryConditions::FullNeumann && 
			Utils::ProgramArgs.Problem.Equation == EquationType::Diffusion && 
			!Utils::ProgramArgs.Problem.ComputeNormalDerivative)
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

		parallelLoop.AggregateChunkResults([this, actions, 
											&A_T_T_Coeffs, &A_T_F_Coeffs, &A_F_F_Coeffs,
											&consistencyCoeffs, &stabilizationCoeffs, &reconstructionCoeffs](AssemblyResult& chunkResult)
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

		this->A_dF_dF = A_F_F.bottomRightCorner(HHO->nDirichletCoeffs, HHO->nDirichletCoeffs);

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
			SparseMatrix tmp = SparseMatrix(A_T_ndF.transpose()) * CellSpace.Solve_A_T_T(A_T_ndF);
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
			Utils::Empty(this->A_dF_dF);
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

	double a(const Vector& u_T, const Vector& u_ndF, const Vector& u_dF, const Vector& v_T, const Vector& v_ndF, const Vector& v_dF)
	{
		Vector Au_T   = A_T_T               * u_T + A_T_ndF              * u_ndF + A_T_dF   * u_dF;
		Vector Au_ndF = A_T_ndF.transpose() * u_T + A_ndF_ndF            * u_ndF + A_ndF_dF * u_dF;
		Vector Au_dF  = A_T_dF.transpose()  * u_T + A_ndF_dF.transpose() * u_ndF + A_dF_dF  * u_dF;
		return Au_T.dot(v_T) + Au_ndF.dot(v_ndF) + Au_dF.dot(v_dF);
	}

	void InitReferenceShapes()
	{
		InitReferenceShapes(this->HHO, &this->TestCase->DiffField);
	}

	// Compute some useful integrals on reference element and store them.
	// Defined for each Dim at the end of the file.
	static void InitReferenceShapes(HHOParameters<Dim>* hho, DiffusionField<Dim>* diffField) { assert(false); }

	//------------------------------//
	//           Init HHO           //
	//------------------------------//
	void InitHHO(bool assembleLocalMatrices = true)
	{
		InitHHO_Faces();
		InitHHO_Elements(assembleLocalMatrices);

		CellSpace             = HHOCellSpace            (_mesh, HHO, _hhoElements);
		ReconstructSpace      = HHOReconstructSpace     (_mesh, HHO, _hhoElements);
		SkeletonSpace         = HHOSkeletonSpace        (_mesh, HHO, _hhoFaces);
		BoundarySpace         = HHOBoundarySpace        (_mesh, HHO, _hhoFaces);
		DirichletSpace        = HHODirichletSpace       (_mesh, HHO, _hhoFaces);
		NeumannSpace          = HHONeumannSpace         (_mesh, HHO, _hhoFaces);
		NonDirichletFaceSpace = HHONonDirichletFaceSpace(_mesh, HHO, _hhoFaces);
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
		return CellSpace.InnerProdWithBasis(sourceFunction);
	}

	Vector AssembleSourceTerm(const Vector& sourceFuncCoeffs)
	{
		if (sourceFuncCoeffs.rows() == HHO->nTotalCellUnknowns) // degree k
			return CellSpace.ApplyMassMatrix(sourceFuncCoeffs);

		if (sourceFuncCoeffs.rows() == HHO->nTotalReconstructUnknowns) // degree k+1
		{
			Vector b_source = Vector(HHO->nTotalCellUnknowns);
			ElementParallelLoop<Dim> parallelLoop(this->_mesh->Elements);
			parallelLoop.Execute([this, &sourceFuncCoeffs, &b_source](Element<Dim>* e)
				{
					Diff_HHOElement<Dim>* element = this->HHOElement(e);
					b_source.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns) = element->ApplyCellReconstructMassMatrix(sourceFuncCoeffs.segment(e->Number * HHO->nReconstructUnknowns, HHO->nReconstructUnknowns));
				}
			);
			return b_source;
		}
		Utils::FatalError("the argument sourceFuncCoeffs does not have a correct size");
		return Vector();
	}

	// Solution on the Dirichlet faces
	Vector AssembleDirichletTerm(DomFunction dirichletFunction)
	{
		return DirichletSpace.Project(dirichletFunction);
	}


	//         [    0    ]   <---- interior faces
	// returns [---------]
	//         [ b_neumF ]   <---- Neumann faces
	// 
	// where b_neumF = [ (neumannFunc|phi_i)_F ]
	Vector AssembleNeumannTerm(DomFunction neumannFunction)
	{
		Vector b_neumann = Vector::Zero(HHO->nTotalFaceUnknowns);
		b_neumann.tail(HHO->nNeumannUnknowns) = NeumannSpace.InnerProdWithBasis(neumannFunction);
		return b_neumann;
	}

	//         [    0    ]   <---- interior faces
	// returns [---------]
	//         [ b_neumF ]   <---- Neumann faces
	// 
	// where b_neumF = [ M * coeffs ]
	Vector AssembleNeumannTerm(const Vector& neumannFuncCoeffs)
	{
		Vector b_neumann = Vector::Zero(HHO->nTotalFaceUnknowns);
		b_neumann.tail(HHO->nNeumannUnknowns) = NeumannSpace.ApplyMassMatrix(neumannFuncCoeffs);
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

	/*Vector SolveCellUnknownsOnBoundaryOnly(const Vector& faceUnknowns, const Vector& b_T)
	{
		BigNumber nBdryCellUnknowns = _mesh->NBoundaryElements() * HHO->nCellUnknowns;
		Vector v = b_T.head(nBdryCellUnknowns) - A_T_ndF.topRows(nBdryCellUnknowns) * faceUnknowns;

		Vector cellUnknowns(nBdryCellUnknowns);

		ElementParallelLoop<Dim> parallelLoop(this->_mesh->BoundaryElements);
		parallelLoop.Execute([this, &v, &cellUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
					Diff_HHOElement<Dim>* elem = HHOElement(e);
					cellUnknowns.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns) = elem->AttSolver.solve(v.segment(e->Number * HHO->nCellUnknowns, HHO->nCellUnknowns));
			});
		return cellUnknowns;
	}*/

	Vector SolveCellUnknowns(const Vector& faceUnknowns, const Vector& b_T)
	{
		assert(faceUnknowns.rows() == HHO->nTotalFaceUnknowns);
		assert(b_T.rows() == HHO->nTotalCellUnknowns);

		return CellSpace.Solve_A_T_T(b_T - A_T_ndF * faceUnknowns);
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
		this->b = this->B_ndF -this->A_T_ndF.transpose() * CellSpace.Solve_A_T_T(this->B_T);
		return this->b;
	}

	Vector CondensedRHS(const Vector& b_T, const Vector& b_ndF)
	{
		assert(b_T.rows() == HHO->nTotalCellUnknowns);
		assert(b_ndF.rows() == HHO->nTotalFaceUnknowns);

		return b_ndF - this->A_T_ndF.transpose() * CellSpace.Solve_A_T_T(b_T);
	}

	Vector CondensedRHS_noNeumannZeroDirichlet(const Vector& b_T)
	{
		assert(b_T.rows() == HHO->nTotalCellUnknowns);

		return -this->A_T_ndF.transpose() * CellSpace.Solve_A_T_T(b_T);
	}

	Vector CondensedRHS_noDirichletZeroNeumann(const Vector& b_T)
	{
		assert(b_T.rows() == HHO->nTotalCellUnknowns);

		return -this->A_T_ndF.transpose() * CellSpace.Solve_A_T_T(b_T);
	}

	//------------------------//
	// For biharmonic problem //
	//------------------------//

	SparseMatrix Theta_T_bF_transpose()
	{
		FaceParallelLoop<Dim> parallelLoop(_mesh->BoundaryFaces);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nCellUnknowns * HHO->nFaceUnknowns);

		parallelLoop.Execute([this](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				Element<Dim>* e = f->Element1;
				int i = f->Number - HHO->nInteriorFaces;
				int j = _mesh->BoundaryElementNumber(e);

				DenseMatrix S = HHOElement(e)->SolveCellUnknownsMatrix().middleCols(e->LocalNumberOf(f) * HHO->nFaceUnknowns, HHO->nFaceUnknowns);

				chunk->Results.Coeffs.Add(i * HHO->nFaceUnknowns, j * HHO->nCellUnknowns, S.transpose());
			});

		SparseMatrix mat(HHO->nBoundaryFaces * HHO->nFaceUnknowns, _mesh->NBoundaryElements() * HHO->nCellUnknowns);
		parallelLoop.Fill(mat);
		return mat;
	}
	
	// Mass matrix of degree k, on boundary elements only
	/*SparseMatrix CellMassMatrixOnBoundaryElements()
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
	}*/

	Vector ExtractElemBoundary(const Vector v)
	{
		assert(_mesh->BoundaryElementsNumberedFirst());
		if (v.rows() == HHO->nTotalReconstructUnknowns)
			return v.head(_mesh->BoundaryElements.size() * HHO->nReconstructUnknowns);
		else if (v.rows() == HHO->nTotalCellUnknowns)
			return v.head(_mesh->BoundaryElements.size() * HHO->nCellUnknowns);
		Utils::FatalError("Not implemented");
		return Vector();
	}

	/*SparseMatrix A_bT_bT_Matrix()
	{
		assert(_mesh->BoundaryElementsNumberedFirst());

		ElementParallelLoop<Dim> parallelLoop(_mesh->BoundaryElements);
		parallelLoop.ReserveChunkCoeffsSize(HHO->nCellUnknowns * HHO->nCellUnknowns);
		parallelLoop.Execute([this](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* elem = HHOElement(e);
				chunk->Results.Coeffs.Add(e->Number * HHO->nCellUnknowns, e->Number * HHO->nCellUnknowns, elem->A.topLeftCorner(HHO->nCellUnknowns, HHO->nCellUnknowns));
			});
		SparseMatrix mat = SparseMatrix(_mesh->BoundaryElements.size() * HHO->nCellUnknowns, _mesh->BoundaryElements.size() * HHO->nCellUnknowns);
		parallelLoop.Fill(mat);
		return mat;
	}*/


	//---------------------------------------//
	//                Exports                //
	//---------------------------------------//

public:
	void ExportReconstructedVectorToGMSH(const Vector& reconstructedVector, const ExportModule& out, string name, double tolerance = 1e-3, int maxRefinements = 6, bool takeAbsoluteValue = false)
	{
		if (HHO->OrthogonalizeElemBases())
			Utils::Error("The export to GMSH has not been implemented when the bases are orthonormalized against each element.");
		this->_mesh->ExportToGMSH_Elements(this->HHO->ReconstructionBasis, reconstructedVector, out.GetFilePathPrefix(), name, tolerance, maxRefinements, takeAbsoluteValue);
	}

	void ExportSolutionToGMSH(const Vector& reconstructedSolution, const ExportModule& out, double tolerance = 1e-3, int maxRefinements = 6)
	{
		ExportReconstructedVectorToGMSH(reconstructedSolution, out, "potential", tolerance, maxRefinements);
	}

	void ExportErrorToGMSH(const Vector& faceCoeffs, const ExportModule& out)
	{
		if (HHO->OrthogonalizeElemBases())
			Utils::Error("The export to GMSH has not been implemented when the bases are orthonormalized against each element.");
		Vector cellCoeffs = CellSpace.Solve_A_T_T(B_T - A_T_ndF * faceCoeffs);
		this->_mesh->ExportToGMSH_Elements(this->HHO->CellBasis, cellCoeffs, out.GetFilePathPrefix(), "error");
	}

private:
	BigNumber FirstDOFGlobalNumber(Diff_HHOElement<Dim>* element)
	{
		return element->Number() * HHO->CellBasis->Size();
	}
	BigNumber FirstDOFGlobalNumber(Diff_HHOFace<Dim>* face)
	{
		return this->HHO->nTotalCellUnknowns + face->Number() * HHO->FaceBasis->Size();
	}
};

#ifdef ENABLE_1D
template <>
void Diffusion_HHO<1>::InitReferenceShapes(HHOParameters<1>* hho, DiffusionField<1>* diffField)
{
	FunctionalBasis<1>* reconstructionBasis = hho->ReconstructionBasis;
	FunctionalBasis<1>* cellBasis = hho->CellBasis;
	FunctionalBasis<0>* faceBasis = hho->FaceBasis;

	// Elements
	if (hho->OrthogonalizeElemBases())
	{
		OrthogonalBasis<1>* orthoBasis = Segment::InitReferenceShape()->Orthogonalize(reconstructionBasis, hho->NElemOrthogonalizations(), hho->OrthonormalizeElemBases());
		Segment::InitReferenceShape()->ComputeAndStoreStiffnessMatrices(orthoBasis);
		Segment::InitReferenceShape()->ComputeAndStoreIntegralVector(orthoBasis);
	}
	else
	{
		// Store mass matrices
		Segment::InitReferenceShape()->ComputeAndStoreMassMatrix(cellBasis);
		Segment::InitReferenceShape()->ComputeAndStoreMassMatrix(reconstructionBasis);
		Segment::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix(cellBasis, reconstructionBasis);
		// Store stiffness matrices
		Segment::InitReferenceShape()->ComputeAndStoreStiffnessMatrices(reconstructionBasis);
		// Integrals
		Segment::InitReferenceShape()->ComputeAndStoreIntegralVector(reconstructionBasis);
	}

	// Faces
	if (hho->OrthogonalizeFaceBases())
		CartesianShape<1, 0>::InitReferenceShape()->Orthogonalize(faceBasis, hho->NFaceOrthogonalizations(), hho->OrthonormalizeFaceBases());
	else
		CartesianShape<1, 0>::InitReferenceShape()->ComputeAndStoreMassMatrix(faceBasis);
}
#endif // ENABLE_1D

#ifdef ENABLE_2D
template <>
void Diffusion_HHO<2>::InitReferenceShapes(HHOParameters<2>* hho, DiffusionField<2>* diffField)
{
	FunctionalBasis<2>* reconstructionBasis = hho->ReconstructionBasis;
	FunctionalBasis<2>* cellBasis = hho->CellBasis;
	FunctionalBasis<1>* faceBasis = hho->FaceBasis;

	// Elements
	if (hho->OrthogonalizeElemBases())
	{
		// No need to store mass matrices.
		// 
		// Store stiffness matrices
		OrthogonalBasis<2>* orthoBasis;
		orthoBasis = Triangle::InitReferenceShape()->Orthogonalize(reconstructionBasis, hho->NElemOrthogonalizations(), hho->OrthonormalizeElemBases());
		Triangle::InitReferenceShape()->ComputeAndStoreStiffnessMatrices(orthoBasis);
		Triangle::InitReferenceShape()->ComputeAndStoreIntegralVector(orthoBasis);

		orthoBasis = Quadrilateral::InitReferenceShape()->Orthogonalize(reconstructionBasis, hho->NElemOrthogonalizations(), hho->OrthonormalizeElemBases());
		Quadrilateral::InitReferenceShape()->ComputeAndStoreStiffnessMatrices(orthoBasis);
		Quadrilateral::InitReferenceShape()->ComputeAndStoreIntegralVector(orthoBasis);

		orthoBasis = CartesianShape<2>::InitReferenceShape()->Orthogonalize(reconstructionBasis, hho->NElemOrthogonalizations(), hho->OrthonormalizeElemBases());
		CartesianShape<2>::InitReferenceShape()->ComputeAndStoreStiffnessMatrices(orthoBasis);
		CartesianShape<2>::InitReferenceShape()->ComputeAndStoreIntegralVector(orthoBasis);
	}
	else
	{
		// Store mass matrices
		Triangle::InitReferenceShape()->ComputeAndStoreMassMatrix(cellBasis);
		Triangle::InitReferenceShape()->ComputeAndStoreMassMatrix(reconstructionBasis);
		Triangle::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix(cellBasis, reconstructionBasis);
		Quadrilateral::InitReferenceShape()->ComputeAndStoreMassMatrix(cellBasis);
		Quadrilateral::InitReferenceShape()->ComputeAndStoreMassMatrix(reconstructionBasis);
		Quadrilateral::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix(cellBasis, reconstructionBasis);
		CartesianShape<2>::InitReferenceShape()->ComputeAndStoreMassMatrix(cellBasis);
		CartesianShape<2>::InitReferenceShape()->ComputeAndStoreMassMatrix(reconstructionBasis);
		CartesianShape<2>::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix(cellBasis, reconstructionBasis);
		// Stiffness matrices
		Triangle::InitReferenceShape()->ComputeAndStoreStiffnessMatrices(reconstructionBasis);
		Quadrilateral::InitReferenceShape()->ComputeAndStoreStiffnessMatrices(reconstructionBasis);
		CartesianShape<2>::InitReferenceShape()->ComputeAndStoreStiffnessMatrices(reconstructionBasis);
		if (diffField)
		{
			for (auto K : diffField->Tensors())
			{
				Quadrilateral::InitReferenceShape()->ComputeAndStoreReconstructStiffnessMatrix(*K, reconstructionBasis);
				CartesianShape<2>::InitReferenceShape()->ComputeAndStoreReconstructStiffnessMatrix(*K, reconstructionBasis);
			}
		}
		// Integrals
		Triangle::InitReferenceShape()->ComputeAndStoreIntegralVector(reconstructionBasis);
		Quadrilateral::InitReferenceShape()->ComputeAndStoreIntegralVector(reconstructionBasis);
		CartesianShape<2>::InitReferenceShape()->ComputeAndStoreIntegralVector(reconstructionBasis);
	}

	// Faces
	if (hho->OrthogonalizeFaceBases())
	{
		Segment::InitReferenceShape()->Orthogonalize(faceBasis, hho->NFaceOrthogonalizations(), hho->OrthonormalizeFaceBases());
		CartesianShape<2, 1>::InitReferenceShape()->Orthogonalize(faceBasis, hho->NFaceOrthogonalizations(), hho->OrthonormalizeFaceBases());
	}
	else
	{
		// Mass matrices
		Segment::InitReferenceShape()->ComputeAndStoreMassMatrix(faceBasis);
		CartesianShape<2, 1>::InitReferenceShape()->ComputeAndStoreMassMatrix(faceBasis);
	}
}
#endif // ENABLE_2D

#ifdef ENABLE_3D
template <>
void Diffusion_HHO<3>::InitReferenceShapes(HHOParameters<3>* hho, DiffusionField<3>* diffField)
{
	FunctionalBasis<3>* reconstructionBasis = hho->ReconstructionBasis;
	FunctionalBasis<3>* cellBasis = hho->CellBasis;
	FunctionalBasis<2>* faceBasis = hho->FaceBasis;

	// Elements
	if (hho->OrthogonalizeElemBases())
	{
		// No need to store mass matrices.
		// 
		// Store stiffness matrices
		OrthogonalBasis<3>* orthoBasis;
		orthoBasis = Tetrahedron::InitReferenceShape()->Orthogonalize(reconstructionBasis, hho->NElemOrthogonalizations(), hho->OrthonormalizeElemBases());
		Tetrahedron::InitReferenceShape()->ComputeAndStoreStiffnessMatrices(orthoBasis);
		Tetrahedron::InitReferenceShape()->ComputeAndStoreIntegralVector(orthoBasis);

		orthoBasis = CartesianShape<3>::InitReferenceShape()->Orthogonalize(reconstructionBasis, hho->NElemOrthogonalizations(), hho->OrthonormalizeElemBases());
		CartesianShape<3>::InitReferenceShape()->ComputeAndStoreStiffnessMatrices(orthoBasis);
		CartesianShape<3>::InitReferenceShape()->ComputeAndStoreIntegralVector(orthoBasis);
	}
	else
	{
		// Store mass matrices
		Tetrahedron::InitReferenceShape()->ComputeAndStoreMassMatrix(cellBasis);
		Tetrahedron::InitReferenceShape()->ComputeAndStoreMassMatrix(reconstructionBasis);
		Tetrahedron::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix(cellBasis, reconstructionBasis);
		CartesianShape<3>::InitReferenceShape()->ComputeAndStoreMassMatrix(cellBasis);
		CartesianShape<3>::InitReferenceShape()->ComputeAndStoreMassMatrix(reconstructionBasis);
		CartesianShape<3>::InitReferenceShape()->ComputeAndStoreCellReconstructMassMatrix(cellBasis, reconstructionBasis);
		// Stiffness matrices
		Tetrahedron::InitReferenceShape()->ComputeAndStoreStiffnessMatrices(reconstructionBasis);
		CartesianShape<3>::InitReferenceShape()->ComputeAndStoreStiffnessMatrices(reconstructionBasis);
		if (diffField)
		{
			for (auto K : diffField->Tensors())
				CartesianShape<3>::InitReferenceShape()->ComputeAndStoreReconstructStiffnessMatrix(*K, reconstructionBasis);
		}
		// Integrals
		Tetrahedron::InitReferenceShape()->ComputeAndStoreIntegralVector(reconstructionBasis);
		CartesianShape<3>::InitReferenceShape()->ComputeAndStoreIntegralVector(reconstructionBasis);
	}

	// Faces
	if (hho->OrthogonalizeFaceBases())
	{
		Triangle::InitReferenceShape()->Orthogonalize(faceBasis, hho->NFaceOrthogonalizations(), hho->OrthonormalizeFaceBases());
		Quadrilateral::InitReferenceShape()->Orthogonalize(faceBasis, hho->NFaceOrthogonalizations(), hho->OrthonormalizeFaceBases());
		CartesianShape<3, 2>::InitReferenceShape()->Orthogonalize(faceBasis, hho->NFaceOrthogonalizations(), hho->OrthonormalizeFaceBases());
	}
	else
	{
		// Mass matrices
		Triangle::InitReferenceShape()->ComputeAndStoreMassMatrix(faceBasis);
		Quadrilateral::InitReferenceShape()->ComputeAndStoreMassMatrix(faceBasis);
		CartesianShape<3, 2>::InitReferenceShape()->ComputeAndStoreMassMatrix(faceBasis);
	}
}
#endif // ENABLE_3D