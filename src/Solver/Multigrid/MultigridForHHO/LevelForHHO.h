#pragma once
#include "../Level.h"
#include "../../../HHO/Diffusion_HHO.h"
#include "../../../HHO/ZeroMeanEnforcer.h"
using namespace std;

template <int Dim>
class LevelForHHO : public Level
{
private:
	GMG_H_Prolongation _hProlongation = GMG_H_Prolongation::CellInterp_Trace;
	GMG_P_Prolongation _pProlongation = GMG_P_Prolongation::Injection;
	GMG_P_Restriction _pRestriction = GMG_P_Restriction::RemoveHigherOrders;
	bool _useHigherOrderReconstruction;
	bool _useHeterogeneousWeighting;
	ZeroMeanEnforcerFromFaceCoeffs<Dim> _zeroMean;
public:
	Diffusion_HHO<Dim>* _problem;

	LevelForHHO(int number, Diffusion_HHO<Dim>* problem,
		GMG_H_Prolongation hProlongation, GMG_P_Prolongation pProlongation, GMG_P_Restriction pRestriction,
		bool useHigherOrderReconstruction, bool useHeterogeneousWeighting)
		:
		Level(number),
		_zeroMean(problem)
	{
		this->_problem = problem;
		this->_hProlongation = hProlongation;
		this->_pProlongation = pProlongation;
		this->_pRestriction = pRestriction;
		this->_useHigherOrderReconstruction = useHigherOrderReconstruction;
		this->_useHeterogeneousWeighting = useHeterogeneousWeighting;
	}

	BigNumber NUnknowns() override
	{
		return _problem->HHO->nTotalFaceUnknowns;
	}

	int PolynomialDegree() override
	{
		return _problem->HHO->FaceBasis->GetDegree();
	}

	int BlockSizeForBlockSmoothers() override
	{
		return _problem->HHO->nFaceUnknowns;
	}

	void CoarsenMesh(H_CoarsStgy elemCoarseningStgy, FaceCoarseningStrategy faceCoarseningStgy, double coarseningFactor, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) override
	{
		Mesh<Dim>* mesh = _problem->_mesh;
		if (Utils::IsRefinementStrategy(elemCoarseningStgy) && mesh->CoarseMesh == nullptr)
		{
			noCoarserMeshProvided = true;
			return;
		}
		cout << "\tCoarsening mesh " << mesh->Id << endl;
		mesh->CoarsenMesh(elemCoarseningStgy, faceCoarseningStgy, coarseningFactor);
		if (!mesh->CoarseMesh || mesh->CoarseMesh->InteriorFaces.size() == 0)
			coarsestPossibleMeshReached = true;
	}

	void ExportMeshToMatlab(Mesh<Dim>* levelMesh)
	{
		string filePath = this->Out.GetFilePath("mesh", ".m");
		levelMesh->ExportToMatlab2(filePath);
	}

private:
	void SetupDiscretizedOperator() override 
	{
		if (this->ComesFrom == CoarseningType::P && _problem->HHO->FaceBasis->IsHierarchical)
		{
			LevelForHHO<Dim>* finerLevel = dynamic_cast<LevelForHHO<Dim>*>(FinerLevel);

			int nHigherDegreeUnknowns = finerLevel->_problem->HHO->nFaceUnknowns;
			int nLowerDegreeUnknowns = this->_problem->HHO->nFaceUnknowns;
			NumberParallelLoop<CoeffsChunk> parallelLoop(finerLevel->OperatorMatrix->rows() / nHigherDegreeUnknowns);
			parallelLoop.Execute([this, nHigherDegreeUnknowns, nLowerDegreeUnknowns](BigNumber i, ParallelChunk<CoeffsChunk>* chunk)
				{
					for (int k = 0; k < nLowerDegreeUnknowns; k++)
					{
						for (RowMajorSparseMatrix::InnerIterator it(*this->FinerLevel->OperatorMatrix, i*nHigherDegreeUnknowns + k); it; ++it)
						{
							auto j = it.col() / nHigherDegreeUnknowns;
							int l = it.col() - j * nHigherDegreeUnknowns;
							if (l < nLowerDegreeUnknowns)
								chunk->Results.Coeffs.Add(i*nLowerDegreeUnknowns + k, j*nLowerDegreeUnknowns + l, it.value());
						}
					}
				});
			SparseMatrix* Ac = new SparseMatrix(this->_problem->HHO->nTotalFaceUnknowns, this->_problem->HHO->nTotalFaceUnknowns);
			parallelLoop.Fill(*Ac);
			this->OperatorMatrix = Ac;
		}
		else
		{
			assert(this->_problem->A.rows() > 0 && "The matrix has not been assembled");
			this->OperatorMatrix = &this->_problem->A;
		}
	}

	/*Smoother* CreateSmoother(string smootherCode, int nSmootherIterations, int blockSize, double omega) override
	{
		if (smootherCode.compare(MegaSmootherInternalSolver<Dim>::Code()) == 0)
			return new MegaSmoother<Dim>(this->_problem->_mesh, blockSize, nSmootherIterations);
		else
			return Level::CreateSmoother(smootherCode, nSmootherIterations, blockSize, omega);
	}*/

	void OnStartSetup() override
	{
		cout << "\t\tk = " << this->PolynomialDegree() << endl;
		if (_hProlongation == GMG_H_Prolongation::FaceInject)
			cout << "\t\tMesh                : " << this->_problem->_mesh->Faces.size() << " faces" << endl;
		else
		{
			cout << "\t\tMesh " << this->_problem->_mesh->Id << "              : " << this->_problem->_mesh->Elements.size() << " elements, regularity = ";
			double regularity = this->_problem->_mesh->Regularity();
			if (regularity == 0)
				cout << "unknown";
			else
				cout << regularity;
			if (!this->IsFinestLevel())
				cout << ", coarsening factor = " << this->_problem->_mesh->CoarseningFactor();
			cout << endl;
		}

		if (ExportComponents)
			this->ExportMeshToMatlab(this->_problem->_mesh);

		if (!IsCoarsestLevel())
		{
			if (this->ComesFrom == CoarseningType::P && 
				_pProlongation == GMG_P_Prolongation::Injection && _pRestriction == GMG_P_Restriction::RemoveHigherOrders && // so the problem hasn't been assembled at this level
				(this->CoarserLevel->ComesFrom == CoarseningType::H || this->CoarserLevel->ComesFrom == CoarseningType::HP)) 
			{
				this->_problem->InitReferenceShapes();
				if (_hProlongation == GMG_H_Prolongation::CellInterp_FinerApproxL2proj_Trace)
					this->_problem->InitHHO(false);
				else
					this->_problem->InitHHO_Faces(); // for the inverses of the face mass matrices (to compute the trace on the faces at the end of the prolongation)
			}

			Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

			if (this->CoarserLevel->ComesFrom == CoarseningType::P && _pProlongation == GMG_P_Prolongation::H_Prolongation)
			{
				coarsePb->InitReferenceShapes();
				coarsePb->InitHHO();
			}
			else if (this->CoarserLevel->ComesFrom == CoarseningType::H || this->CoarserLevel->ComesFrom == CoarseningType::HP)
			{
				if (this->UseGalerkinOperator)
					coarsePb->InitHHO();
				else
				{
					ActionsArguments actions;
					actions.LogAssembly = false;
					actions.AssembleRightHandSide = false;
					//actions.InitReferenceShapes = (CoarserLevel->ComesFrom == CoarseningType::HP || CoarserLevel->ComesFrom == CoarseningType::H)
						//&& this->ComesFrom == CoarseningType::P;
					actions.InitReferenceShapes = true;
					coarsePb->Assemble(actions);
				}
			}
		}
	}

	void SetupProlongation() override
	{
		if (this->CoarserLevel->ComesFrom == CoarseningType::H || this->CoarserLevel->ComesFrom == CoarseningType::HP)
			SetupProlongation_H();
		else if (this->CoarserLevel->ComesFrom == CoarseningType::P)
			SetupProlongation_P();
		else
			assert(false);
	}

	void SetupProlongation_H()
	{
		Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;


		/*cout << endl;
		Timer totalProlongTimer;
		totalProlongTimer.Start();*/

		if (_hProlongation == GMG_H_Prolongation::CellInterp_Trace)
		{
			//----------------------------------------------------------//
			//                       Nested meshes                      //
			// Step 1: Interpolation from coarse faces to coarse cells. //
			// Step 2: Trace on the fine faces.                         //
			//----------------------------------------------------------//

			SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
			SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCoarseCellsToFineFaces();

			if (ExportComponents)
			{
				Level::ExportMatrix(I_c, "I_c");
				Level::ExportMatrix(Pi_f, "Pi_f");

				SparseMatrix SolveCellUnknowns = GetSolveCellUnknownsMatrix(coarsePb);
				Level::ExportMatrix(SolveCellUnknowns, "SolveCellUnknowns");
			}

			P = Pi_f * I_c;
		}
		else if (_hProlongation == GMG_H_Prolongation::CellInterp_Inject_Trace)
		{
			//----------------------------------------------------------//
			//          Same method, other implementation               //
			// Step 1: Interpolation from coarse faces to coarse cells. //
			// Step 2: Canonical injection from coarse to fine cells.   //
			// Step 3: Trace on the fine faces.                         //
			//----------------------------------------------------------//

			SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
			SparseMatrix J_f_c = GetGlobalCanonicalInjectionMatrixCoarseToFineElements();
			SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCellsToFaces(finePb);

			if (ExportComponents)
			{
				Level::ExportMatrix(I_c, "I_c");
				Level::ExportMatrix(J_f_c, "J_f_c");
				Level::ExportMatrix(Pi_f, "Pi_f");

				SparseMatrix SolveCellUnknowns = GetSolveCellUnknownsMatrix(coarsePb);
				Level::ExportMatrix(SolveCellUnknowns, "SolveCellUnknowns");
			}

			P = Pi_f * J_f_c * I_c;
		}
		else if (_hProlongation == GMG_H_Prolongation::CellInterp_ExactL2proj_Trace)
		{
			//---------------------------------------------------------------------//
			//                      Non-nested variant                             //
			// Step 1: Interpolation from coarse faces to coarse cells.            //
			// Step 2: Instead of the canonical injection which is now impossible, //
			//         L2-projection onto the fine cells.                          //
			// Step 3: Trace on the fine faces.                                    //
			//                                                                     //
			// Details: To be computed exactly, the L2-projection requires the     //
			//          computation of the intersection between the coarse and     //
			//          fine elements, which is costly.                            //
			//---------------------------------------------------------------------//
			
			//Timer intersectionTimer;
			//intersectionTimer.Start();
			coarsePb->_mesh->SetOverlappingFineElementsViaExactIntersection();
			//intersectionTimer.Stop();
			//cout << "Intersections: " << intersectionTimer.CPU().InMilliseconds << endl;

			SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
			//Timer innerProdTimer;
			//innerProdTimer.Start();
			SparseMatrix L2Proj = GetGlobalL2ProjectionMatrixCoarseToFineElements();
			//innerProdTimer.Stop();
			//cout << "Inner products: " << innerProdTimer.CPU().InMilliseconds << endl;
			SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCellsToFaces(finePb);

			coarsePb->_mesh->DeleteOverlappingFineElementsInformation();

			if (ExportComponents)
			{
				Level::ExportMatrix(I_c, "I_c");
				Level::ExportMatrix(L2Proj, "J_f_c");
				Level::ExportMatrix(Pi_f, "Pi_f");
			}

			P = Pi_f * L2Proj * I_c;
		}
		else if (_hProlongation == GMG_H_Prolongation::CellInterp_ApproxL2proj_Trace)
		{
			//---------------------------------------------------------------------//
			//                        Non-nested variant                           //
			//            with approximate L2-projection (less costly)             //
			//                                                                     //
			// Step 1: Interpolation from coarse faces to coarse cells.            //
			// Step 2: Instead of the canonical injection which is now impossible, //
			//         L2-projection onto the fine cells.                          //
			//         However, because the exact L2-proj is costly to compute,    //
			//         we do it approximately in such a way that it has the same   //
			//         approximation properties.                                   //
			// Step 3: Trace on the fine faces.                                    //
			//                                                                     //
			// Details: This is actually the exact same code as                    //
			//          CellInterp_Inject_Trace (the L2-proj is implemented as the //
			//          canonical injection).                                      //
			//---------------------------------------------------------------------//

			// The coarse element associated to every mesh element must be the closest one.
			// Every fine element must appear in only one coarse element's FinerElements list.


			if (!Utils::BuildsNestedMeshHierarchy(coarsePb->_mesh->ComesFrom.CS) && coarsePb->_mesh->ComesFrom.CS != H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours)
			{
				if (Dim == 2 && this->PolynomialDegree() > 1)
					Utils::Warning("The degree k=" + to_string(this->PolynomialDegree()) + " is too high to ensure enough accuracy of the approximate L2-projection. Non-optimal convergence may be observed. Option '-prolong " + to_string((unsigned)GMG_H_Prolongation::CellInterp_FinerApproxL2proj_Trace) + "' (approx. with subdivision) is recommended.");
				else if (Dim == 3)
					Utils::Warning("This rough approximation of the L2-projection does not work well in 3D. Non-optimal convergence may be observed. Option '-prolong " + to_string((unsigned)GMG_H_Prolongation::CellInterp_FinerApproxL2proj_Trace) + "' (approx. with subdivision) is recommended.");
			}

#ifndef NDEBUG
			//coarsePb->_mesh->PlotClustersForApproxL2(false);
#endif // !NDEBUG

			SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
			//Timer innerProdTimer;
			//innerProdTimer.Start();
			SparseMatrix L2Proj = GetGlobalCanonicalInjectionMatrixCoarseToFineElements();
			//innerProdTimer.Stop();
			//cout << "Inner products: " << innerProdTimer.CPU().InMilliseconds << endl;
			SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCellsToFaces(finePb);

			if (ExportComponents)
			{
				Level::ExportMatrix(I_c, "I_c");
				Level::ExportMatrix(L2Proj, "J_f_c");
				Level::ExportMatrix(Pi_f, "Pi_f");
			}

			P = Pi_f * L2Proj * I_c;
		}
		else if (_hProlongation == GMG_H_Prolongation::CellInterp_FinerApproxL2proj_Trace)
		{
			//---------------------------------------------------------------------//
			//                        Non-nested variant                           //
			//               with a finer approximate L2-projection                //
			//                                                                     //
			// Step 1: Interpolation from coarse faces to coarse cells.            //
			// Step 2: Approximate L2-projection onto the fine cells.              //
			// Step 3: Trace on the fine faces.                                    //
			//                                                                     //
			// Details: Derived from the previous approximate L2-projection.       //
			//          The fine elements are now refined and the subelements      //
			//          associated to coarse elements accordingly.                 //
			//          The implementation is the same as the exact L2-proj        //
			//          though. Except that the coarse/fine intersection are the   //
			//          subelements.                                               //
			//---------------------------------------------------------------------//

			if (!Utils::BuildsNestedMeshHierarchy(coarsePb->_mesh->ComesFrom.CS) && coarsePb->_mesh->ComesFrom.CS != H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours)
			{
				if (((Dim == 2 && this->PolynomialDegree() > 3) || (Dim == 3 && this->PolynomialDegree() > 2)) && Utils::ProgramArgs.Solver.MG.NSubtriangulationsForApproxL2Proj == 1)
					Utils::Warning("The degree k=" + to_string(this->PolynomialDegree()) + " is too high to ensure enough accuracy of the approximate L2-projection. Non-optimal convergence may be observed. Consider more subtriangulation.");
			}

			//Timer subdivisionsTimer;
			//subdivisionsTimer.Start();
			coarsePb->_mesh->SetOverlappingFineElementsSubTriangles();
			//subdivisionsTimer.Stop();
			//cout << "Subdivisions: " << subdivisionsTimer.CPU().InMilliseconds << endl;
#ifndef NDEBUG
			//coarsePb->_mesh->PlotClustersForApproxL2(true);
#endif // !NDEBUG

			SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
			//Timer innerProdTimer;
			//innerProdTimer.Start();
			SparseMatrix L2Proj = GetGlobalL2ProjectionMatrixCoarseToFineElements();
			//innerProdTimer.Stop();
			//cout << "Inner products: " << innerProdTimer.CPU().InMilliseconds << endl;
			SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCellsToFaces(finePb);

			if (ExportComponents)
			{
				Level::ExportMatrix(I_c, "I_c");
				Level::ExportMatrix(L2Proj, "J_f_c");
				Level::ExportMatrix(Pi_f, "Pi_f");
			}

			P = Pi_f * L2Proj * I_c;
		}
		else if (_hProlongation == GMG_H_Prolongation::CellInterp_InjectAndTrace)
		{
			//---------------------------------------------------------------------------------------------//
			// Step 1: Interpolation from coarse faces to coarse cells.                                    //
			// Step 2: On faces present on both fine and coarse meshes, we keep the polynomials identical. //
			//         On faces interior to coarse elements, trace of the cell polynomials.                //
			//---------------------------------------------------------------------------------------------//

			SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
			SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCoarseCellsToFineFaces();

			SparseMatrix J_faces = GetGlobalCanonicalInjectionMatrixCoarseToFineFaces();

			if (ExportComponents)
			{
				Level::ExportMatrix(I_c, "I_c");
				Level::ExportMatrix(Pi_f, "Pi_f");
				Level::ExportMatrix(J_faces, "J_faces");

				SparseMatrix SolveCellUnknowns = GetSolveCellUnknownsMatrix(coarsePb);
				Level::ExportMatrix(SolveCellUnknowns, "SolveCellUnknowns");
			}

			SparseMatrix P_algo1 = Pi_f * I_c;
			
			int nFaceUnknowns = finePb->HHO->nFaceUnknowns;

			FaceParallelLoop<Dim> parallelLoop(finePb->_mesh->Faces);
			parallelLoop.ReserveChunkCoeffsSize(nFaceUnknowns * 2 * nFaceUnknowns);

			parallelLoop.Execute([this, P_algo1, J_faces, nFaceUnknowns](Face<Dim>* fineFace, ParallelChunk<CoeffsChunk>* chunk)
				{
					if (fineFace->HasDirichletBC())
						return;

					if (fineFace->IsRemovedOnCoarserGrid)
						chunk->Results.Coeffs.CopyRows(fineFace->Number*nFaceUnknowns, nFaceUnknowns, P_algo1);
					else
						chunk->Results.Coeffs.CopyRows(fineFace->Number*nFaceUnknowns, nFaceUnknowns, J_faces);
				});

			P = SparseMatrix(finePb->HHO->nInteriorAndNeumannFaces * nFaceUnknowns, coarsePb->HHO->nInteriorAndNeumannFaces * nFaceUnknowns);
			parallelLoop.Fill(P);
		}
		else if (_hProlongation == GMG_H_Prolongation::CellInterp_Inject_Adjoint)
		{
			//-------------------------------------------------------------//
			//                           Daniel's                          //
			// Step 1: Interpolation from coarse faces to coarse cells.    //
			// Step 2: Canonical injection from coarse to fine cells.      //
			// Step 2: Adjoint of the same interpolation on the fine mesh. //
			//-------------------------------------------------------------//

			SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
			SparseMatrix J_f_c = GetGlobalCanonicalInjectionMatrixCoarseToFineElements();
			SparseMatrix K_f = GetGlobalMatrixFindFacesWhichReconstructCells(finePb);

			if (ExportComponents)
			{
				Level::ExportMatrix(I_c, "I_c");
				Level::ExportMatrix(J_f_c, "J_f_c");
				Level::ExportMatrix(K_f, "K_f");
			}

			P = K_f * J_f_c * I_c;
		}
		else if (_hProlongation == GMG_H_Prolongation::Wildey)
		{
			//-------------------------------------------------------------------------------------------------//
			//                                           Wildey et al.                                         //
			// The coarse level is built by static condensation of the fine faces interior to coarse elements. //
			// The prolongation solves those condensed unknowns.                                               //
			//-------------------------------------------------------------------------------------------------//

			int nFaceUnknowns = finePb->HHO->nFaceUnknowns;

			ElementParallelLoop<Dim> parallelLoop(coarsePb->_mesh->Elements);

			parallelLoop.Execute([this, finePb, coarsePb, nFaceUnknowns](Element<Dim>* coarseElem, ParallelChunk<CoeffsChunk>* chunk)
				{
					DenseMatrix resolveCondensedFinerFacesFromCoarseBoundary = StaticallyCondenseInteriorFinerFaces(coarseElem, *this->OperatorMatrix, coarsePb, finePb);

					for (int i = 0; i < coarseElem->FinerFacesRemoved.size(); i++)
					{
						Face<Dim>* condensedFineFace = coarseElem->FinerFacesRemoved[i];
						for (Face<Dim>* coarseFace : coarseElem->Faces)
						{
							if (coarseFace->HasDirichletBC())
								continue;

							BigNumber coarseFaceLocalNumberInCoarseElem = coarseElem->LocalNumberOf(coarseFace);
							chunk->Results.Coeffs.Add(condensedFineFace->Number*nFaceUnknowns, coarseFace->Number*nFaceUnknowns, resolveCondensedFinerFacesFromCoarseBoundary.block(i*nFaceUnknowns, coarseFaceLocalNumberInCoarseElem*nFaceUnknowns, nFaceUnknowns, nFaceUnknowns));
						}
					}

					for (Face<Dim>* coarseFace : coarseElem->Faces)
					{
						if (coarseFace->HasDirichletBC())
							continue;

						DenseMatrix local_J_f_c = this->ComputeCanonicalInjectionMatrixCoarseToFine(coarseFace, coarsePb, finePb);
						for (auto fineFace : coarseFace->FinerFaces)
						{
							BigNumber fineFaceLocalNumberInCoarseFace = coarseFace->LocalNumberOf(fineFace);
							chunk->Results.Coeffs.Add(fineFace->Number*nFaceUnknowns, coarseFace->Number*nFaceUnknowns, local_J_f_c.block(fineFaceLocalNumberInCoarseFace*nFaceUnknowns, 0, nFaceUnknowns, nFaceUnknowns));
						}
					}
				});

			P = SparseMatrix(finePb->HHO->nInteriorAndNeumannFaces * nFaceUnknowns, coarsePb->HHO->nInteriorAndNeumannFaces * nFaceUnknowns);
			parallelLoop.Fill(P);
		}
		else if (_hProlongation == GMG_H_Prolongation::FaceInject)
		{
			//----------------------------------------------------------------//
			// Canonical injection from coarse faces to fine faces.           //
			// Implemented to be used with the face coarsening (option -cs f) //
			//----------------------------------------------------------------//

			SparseMatrix J_faces = GetGlobalCanonicalInjectionMatrixCoarseToFineFaces();
			P = J_faces;
		}
		else
			Utils::FatalError("Unknown h-prolongation.");

		//totalProlongTimer.Stop();
		//cout << "Total prolong: " << totalProlongTimer.CPU().InMilliseconds << endl;
	}


	void SetupProlongation_P()
	{
		if (_pProlongation == GMG_P_Prolongation::Injection)
		{
			if (!this->_problem->HHO->FaceBasis->IsHierarchical || !this->_problem->HHO->OrthogonalizeFaceBases())
				Utils::Warning("The natural injection and restriction for p-multigrid are implemented based on the assumption that the face bases are hierarchical and orthogonalized. Degraded convergence may be experienced.");
			// nothing to do
		}
		else if (_pProlongation == GMG_P_Prolongation::H_Prolongation)
		{
			Diffusion_HHO<Dim>* finePb = this->_problem;
			Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

			if (!_useHigherOrderReconstruction)
				Utils::FatalError("The higher-order reconstruction must be enabled to use this p-prolongation operator.");

			SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
			SparseMatrix Pi_f = GetGlobalProjectorMatrixFromHighOrderCellsToFaces(coarsePb, finePb);
			P = Pi_f * I_c;
		}
		else
			Utils::FatalError("Unknown p-prolongation.");
	}


	void SetupRestriction() override
	{
		if (this->CoarserLevel->ComesFrom == CoarseningType::P && _pRestriction == GMG_P_Restriction::RemoveHigherOrders)
		{
			// nothing to do
		}
		else
			R = (RestrictionScalingFactor() * P.transpose()).eval();
	}

	void OnEndSetup() override
	{
		if (_problem->TestCase->BC.Type == PbBoundaryConditions::FullNeumann)
			ComputeVectorsForOrthogonalityConditions();

		// Delete now useless stuff
		if (this->ComesFrom == CoarseningType::H || this->ComesFrom == CoarseningType::HP)
		{
			// If finest level, delete everything you don't need to reconstruct the solution at the end
			if (IsFinestLevel())
			{
				ElementParallelLoop<Dim> parallelLoopE(_problem->_mesh->Elements);
				parallelLoopE.Execute([this](Element<Dim>* element)
					{
						Diff_HHOElement<Dim>* hhoElement = _problem->HHOElement(element);
						hhoElement->DeleteUselessMatricesAfterMultigridSetup();
					});

				FaceParallelLoop<Dim> parallelLoopF(_problem->_mesh->Faces);
				parallelLoopF.Execute([this](Face<Dim>* face)
					{
						Diff_HHOFace<Dim>* hhoFace = _problem->HHOFace(face);
						hhoFace->DeleteUselessMatricesAfterMultigridSetup();
					});
			}

			// If _problem->_mesh is not the fine mesh, delete it if possible
			if (_problem->_mesh->FineMesh && _problem->_mesh != _problem->_mesh->CoarseMesh && Utils::ProgramArgs.Solver.MG.HP_CS != HP_CoarsStgy::H_then_P && Utils::ProgramArgs.Solver.MG.HP_CS != HP_CoarsStgy::HP_then_P)
			{
				_problem->_mesh->FineMesh->CoarseMesh = nullptr;
				_problem->_mesh->CoarseMesh = nullptr;
				_problem->DeleteHHOElements();
				_problem->DeleteHHOFaces();

				//cout << "\tDeleting mesh " << _problem->_mesh->Id << endl;
				//delete _problem->_mesh; // !!!!!!! This delete causes a bug later in the process for large problems, but I can't find out why
			}
		}
	}

public:
	Vector Prolong(Vector& vectorOnTheCoarserLevel) override
	{
		if (this->CoarserLevel->ComesFrom == CoarseningType::P && _pProlongation == GMG_P_Prolongation::Injection)
		{
			// This is possible only if the basis is hierarchical
			Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
			auto nFaces = _problem->HHO->nInteriorAndNeumannFaces;
			auto nHigherDegreeUnknowns = _problem->HHO->nFaceUnknowns;
			auto nLowerDegreeUnknowns = coarsePb->HHO->nFaceUnknowns;
			Vector vectorOnThisLevel = Vector::Zero(nFaces * nHigherDegreeUnknowns);
			for (BigNumber i = 0; i < nFaces; i++)
				vectorOnThisLevel.segment(i*nHigherDegreeUnknowns, nLowerDegreeUnknowns) = vectorOnTheCoarserLevel.segment(i*nLowerDegreeUnknowns, nLowerDegreeUnknowns);
			return vectorOnThisLevel;
		}
		else
			return Level::Prolong(vectorOnTheCoarserLevel);
	}

	Vector Restrict(Vector& vectorOnThisLevel) override
	{
		if (this->CoarserLevel->ComesFrom == CoarseningType::P && _pRestriction == GMG_P_Restriction::RemoveHigherOrders)
		{
			// This makes sense only if the basis is hierarchical and orthogonal
			Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
			auto nFaces = _problem->HHO->nInteriorAndNeumannFaces;
			auto nHigherDegreeUnknowns = _problem->HHO->nFaceUnknowns;
			auto nLowerDegreeUnknowns = coarsePb->HHO->nFaceUnknowns;
			Vector vectorOnTheCoarserLevel(nFaces * nLowerDegreeUnknowns);
			for (BigNumber i = 0; i < nFaces; i++)
				vectorOnTheCoarserLevel.segment(i*nLowerDegreeUnknowns, nLowerDegreeUnknowns) = vectorOnThisLevel.segment(i*nHigherDegreeUnknowns, nLowerDegreeUnknowns);
			return vectorOnTheCoarserLevel;
		}
		else
			return Level::Restrict(vectorOnThisLevel);
	}

	Flops ProlongCost() override
	{
		if (this->CoarserLevel->ComesFrom == CoarseningType::P && _pProlongation == GMG_P_Prolongation::Injection)
			return 0;
		else
			return Level::ProlongCost();
	}

	Flops RestrictCost() override
	{
		if (this->CoarserLevel->ComesFrom == CoarseningType::P && _pRestriction == GMG_P_Restriction::RemoveHigherOrders)
			return 0;
		else
			return Level::RestrictCost();
	}

	//----------------------------------//
	// Full Neumann boundary conditions //
	//----------------------------------//

	// --- Zero-mean condition (to enforce unicity of the solution)
	// 
	// --- Compatibility condition (to ensure existence of a solution)
	//
	// Considering Ax=b with A singular, a solution exists iff b \in Im(A).
	// Fundamental theorem of linear algebra: Im(A) is the orthogonal of Ker(A^T).
	// As A is symmetric, Ker(A^T) = Ker(A), so Im(A) is the orthogonal of Ker(A).
	// Ker(A) is of dimension 1, spanned by the constant functions,
	// so b \in Im(A) iff b orthogonal to the constant function 1.
	// Consequently, if 

private:
	void ComputeVectorsForOrthogonalityConditions()
	{
		_zeroMean.Setup();
	}

public:
	// Fixes the constant in order to have a well-posed problem.
	// Implements (eq.2) above.
	void ApplyZeroMeanCondition(Vector& x) override
	{
		_zeroMean.Enforce(x);
	}
	Flops ApplyZeroMeanConditionCost(Vector& x) override
	{
		return Cost::Dot(x) + Cost::VectorDAXPY(x);
	}


	// Compatibility condition (to enforce the existence of a solution)
	void EnforceCompatibilityCondition(Vector& b) override
	{
		//assert(_gamma.rows() > 0 && _one.rows() > 0 && "ComputeVectorsForOrthogonalityConditions() has not been called!");
		//b -= b.dot(_gamma) * _one;
		_zeroMean.ProjectOntoImage(b);
		assert(_zeroMean.Check(b));
	}
	Flops EnforceCompatibilityConditionCost(Vector& b) override
	{
		//return Cost::Dot(b) + Cost::VectorDAXPY(b);
		return 0;
	}

private:
	double RestrictionScalingFactor()
	{
		double scalingFactor = 1;

		//------------- Try 1 -------------//
		// scalingFactor = 1 / (average eigenvalue of P^T*P)
		/*SparseMatrix PT_P = P.transpose()*P;
		double averageEigenvalue = (PT_P).diagonal().sum() / PT_P.cols();
		//double maxEigenvalue = (PT_P).diagonal().max()
		cout << "averageEigenvalue = " << averageEigenvalue << endl;
		scalingFactor = 1 / averageEigenvalue;*/

		//------------- Try 2 -------------//
		// scalingFactor = coarseSkeletonMeasure / fineSkeletonMeasure;
		/*Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		scalingFactor = coarsePb->_mesh->SkeletonMeasure() / finePb->_mesh->SkeletonMeasure();*/

		//------------- Try 3 -------------//
		// Scaling conservation:
		// Actually rescales correctly, but that's not what we want //
		/*Diffusion_HHO<Dim>* finePb = this->_problem;
		DomFunction constantOne = [](const DomPoint& p) { return 1; };
		Vector constantOneDoFs = finePb->ProjectOnFaceDiscreteSpace(constantOne);

		double skeletonMeasure = 0;
		for (Face<Dim>* face : finePb->_mesh->Faces)
		{
			if (!face->HasDirichletBC())
				skeletonMeasure += face->Measure();
		}

		Vector restrictAndProlongateConstantOne = P * P.transpose() * constantOneDoFs;
		double integral = finePb->IntegralFromFaceDoFs(restrictAndProlongateConstantOne, 0);

		scalingFactor = skeletonMeasure / integral;*/

		//------------- Try 4 -------------//
		/*Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		DomFunction constantOne = [](const DomPoint& p) { return 1; };
		Vector constantOneDoFs = coarsePb->ProjectOnFaceDiscreteSpace(constantOne);

		double coarseProduct = constantOneDoFs.transpose() * coarsePb->A * constantOneDoFs;
		double fineProduct = constantOneDoFs.transpose() * P.transpose() * finePb->A * P * constantOneDoFs;
		scalingFactor = (coarseProduct / fineProduct);

		cout << "scalingFactor = " << scalingFactor << endl;*/

		//------------- Try 5 -------------//
		// scalingFactor = h_coarse / h_fine;
		/*Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		scalingFactor = coarsePb->_mesh->H() / finePb->_mesh->H();*/

		return scalingFactor;
	}

	inline double Weight(Element<Dim>* element, Face<Dim>* face)
	{
		if (face->IsDomainBoundary)
			return 1;
		else if (_useHeterogeneousWeighting)
		{
			auto n = element->OuterNormalVector(face);
			double k = (element->DiffTensor() * n).dot(n);

			auto n1 = face->Element1->OuterNormalVector(face);
			double k1 = (face->Element1->DiffTensor() * n1).dot(n1);

			auto n2 = face->Element2->OuterNormalVector(face);
			double k2 = (face->Element2->DiffTensor() * n2).dot(n2);

			return k / (k1 + k2);
		}
		else
			return 0.5;
		assert(false);
	}

	SparseMatrix GetGlobalInterpolationMatrixFromFacesToCells(Diffusion_HHO<Dim>* problem)
	{
		int nCellUnknowns = _useHigherOrderReconstruction ? problem->HHO->nReconstructUnknowns : problem->HHO->nCellUnknowns;
		int nFaceUnknowns = problem->HHO->nFaceUnknowns;

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, problem, nFaceUnknowns, nCellUnknowns](Element<Dim>* element, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* hhoElement = problem->HHOElement(element);

				DenseMatrix cellInterpMatrix;
				if (_useHigherOrderReconstruction)
					cellInterpMatrix = hhoElement->ReconstructionFromFacesMatrix();
				else
					cellInterpMatrix = hhoElement->SolveCellUnknownsMatrix();

				for (auto face : element->Faces)
				{
					if (face->HasDirichletBC())
						continue;

					BigNumber elemGlobalNumber = element->Number;
					BigNumber faceGlobalNumber = face->Number;
					BigNumber faceLocalNumber = element->LocalNumberOf(face);

					chunk->Results.Coeffs.Add(elemGlobalNumber * nCellUnknowns, faceGlobalNumber * nFaceUnknowns, cellInterpMatrix.block(0, faceLocalNumber*nFaceUnknowns, nCellUnknowns, nFaceUnknowns));
				}
			});

		SparseMatrix M(problem->HHO->nElements * nCellUnknowns, problem->HHO->nTotalFaceUnknowns);
		parallelLoop.Fill(M);

		return M;
	}

	SparseMatrix GetSolveCellUnknownsMatrix(Diffusion_HHO<Dim>* problem)
	{
		int nCellUnknowns = problem->HHO->nCellUnknowns;
		int nFaceUnknowns = problem->HHO->nFaceUnknowns;

		//FunctionalBasis<Dim - 1>* faceBasis = problem->HHO->FaceBasis;
		//FunctionalBasis<Dim>* cellInterpolationBasis = _problem->HHO->CellBasis;

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, problem, nFaceUnknowns, nCellUnknowns](Element<Dim>* element, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* hhoElement = problem->HHOElement(element);
				DenseMatrix reconstructMatrix = hhoElement->SolveCellUnknownsMatrix();

				for (auto face : element->Faces)
				{
					if (face->HasDirichletBC())
						continue;

					BigNumber elemGlobalNumber = element->Number;
					BigNumber faceGlobalNumber = face->Number;
					BigNumber faceLocalNumber = element->LocalNumberOf(face);

					chunk->Results.Coeffs.Add(elemGlobalNumber * nCellUnknowns, faceGlobalNumber * nFaceUnknowns, reconstructMatrix.block(0, faceLocalNumber*nFaceUnknowns, nCellUnknowns, nFaceUnknowns));
				}
			});

		SparseMatrix M(problem->HHO->nElements * nCellUnknowns, problem->HHO->nTotalFaceUnknowns);
		parallelLoop.Fill(M);

		return M;
	}

	SparseMatrix GetGlobalProjectorMatrixFromCellsToFaces(Diffusion_HHO<Dim>* problem)
	{
		int nCellUnknowns = _useHigherOrderReconstruction ? problem->HHO->nReconstructUnknowns : problem->HHO->nCellUnknowns;
		int nFaceUnknowns = problem->HHO->nFaceUnknowns;

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, problem, nCellUnknowns, nFaceUnknowns](Element<Dim>* element, ParallelChunk<CoeffsChunk>* chunk)
			{
				BigNumber elemGlobalNumber = element->Number;
				Diff_HHOElement<Dim>* hhoElem = problem->HHOElement(element);
				FunctionalBasis<Dim>* cellBasis = _useHigherOrderReconstruction ? hhoElem->ReconstructionBasis : hhoElem->CellBasis;
				for (auto face : element->Faces)
				{
					if (face->HasDirichletBC())
						continue;

					BigNumber faceGlobalNumber = face->Number;
					double weight = Weight(element, face);
					Diff_HHOFace<Dim>* hhoFace = problem->HHOFace(face);
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, elemGlobalNumber*nCellUnknowns, weight*hhoFace->Trace(element, cellBasis));
				}
			});

		SparseMatrix Pi(problem->HHO->nTotalFaceUnknowns, problem->HHO->nElements * nCellUnknowns);
		parallelLoop.Fill(Pi);

		return Pi;
	}

	SparseMatrix GetGlobalProjectorMatrixFromHighOrderCellsToFaces(Diffusion_HHO<Dim>* cellProblem, Diffusion_HHO<Dim>* faceProblem)
	{
		int nCellUnknowns = cellProblem->HHO->nReconstructUnknowns;
		int nFaceUnknowns = faceProblem->HHO->nFaceUnknowns;

		ElementParallelLoop<Dim> parallelLoop(cellProblem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, cellProblem, faceProblem, nCellUnknowns, nFaceUnknowns](Element<Dim>* element, ParallelChunk<CoeffsChunk>* chunk)
			{
				BigNumber elemGlobalNumber = element->Number;
				Diff_HHOElement<Dim>* hhoElem = cellProblem->HHOElement(element);
				FunctionalBasis<Dim>* cellBasis = hhoElem->ReconstructionBasis;
				for (auto face : element->Faces)
				{
					if (face->HasDirichletBC())
						continue;

					BigNumber faceGlobalNumber = face->Number;
					double weight = Weight(element, face);
					Diff_HHOFace<Dim>* hhoFace = faceProblem->HHOFace(face);
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, elemGlobalNumber*nCellUnknowns, weight*hhoFace->Trace(element, cellBasis));
				}
			});

		SparseMatrix Pi(faceProblem->HHO->nTotalFaceUnknowns, cellProblem->HHO->nElements * nCellUnknowns);
		parallelLoop.Fill(Pi);

		return Pi;
	}

	SparseMatrix GetGlobalProjectorMatrixFromCoarseCellsToFineFaces()
	{
		Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

		int nCellUnknowns = _useHigherOrderReconstruction ? coarsePb->HHO->nReconstructUnknowns : coarsePb->HHO->nCellUnknowns;
		int nFaceUnknowns = finePb->HHO->nFaceUnknowns;

		FaceParallelLoop<Dim> parallelLoop(finePb->_mesh->Faces);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, finePb, coarsePb, nCellUnknowns, nFaceUnknowns](Face<Dim>* face, ParallelChunk<CoeffsChunk>* chunk)
			{
				if (face->HasDirichletBC())
					return;

				Diff_HHOFace<Dim>* hhoFace = finePb->HHOFace(face);

				BigNumber faceGlobalNumber = face->Number;
				if (face->IsDomainBoundary || face->IsRemovedOnCoarserGrid)
				{
					Element<Dim>* coarseElem = face->Element1->CoarserElement;
					BigNumber coarseElemGlobalNumber = coarseElem->Number;
					Diff_HHOElement<Dim>* hhoCoarseElem = coarsePb->HHOElement(coarseElem);
					FunctionalBasis<Dim>* cellInterpolationBasis = _useHigherOrderReconstruction ? hhoCoarseElem->ReconstructionBasis : hhoCoarseElem->CellBasis;

					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, coarseElemGlobalNumber*nCellUnknowns, hhoFace->Trace(coarseElem, cellInterpolationBasis));
				}
				else
				{
					Element<Dim>* coarseElem1 = face->Element1->CoarserElement;
					Element<Dim>* coarseElem2 = face->Element2->CoarserElement;
					BigNumber coarseElem1GlobalNumber = coarseElem1->Number;
					BigNumber coarseElem2GlobalNumber = coarseElem2->Number;
					Diff_HHOElement<Dim>* hhoCoarseElem1 = coarsePb->HHOElement(coarseElem1);
					Diff_HHOElement<Dim>* hhoCoarseElem2 = coarsePb->HHOElement(coarseElem2);
					FunctionalBasis<Dim>* cellInterpolationBasis1 = _useHigherOrderReconstruction ? hhoCoarseElem1->ReconstructionBasis : hhoCoarseElem1->CellBasis;
					FunctionalBasis<Dim>* cellInterpolationBasis2 = _useHigherOrderReconstruction ? hhoCoarseElem2->ReconstructionBasis : hhoCoarseElem2->CellBasis;

					//double weight1 = Weight(coarseElem1, face->CoarseFace); // Careful with using face->CoarseFace when the mesh isn't nested...
					double weight1 = Weight(face->Element1, face);
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, coarseElem1GlobalNumber*nCellUnknowns, weight1*hhoFace->Trace(coarseElem1, cellInterpolationBasis1));

					//double weight2 = Weight(coarseElem2, face->CoarseFace);
					double weight2 = Weight(face->Element2, face);
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, coarseElem2GlobalNumber*nCellUnknowns, weight2*hhoFace->Trace(coarseElem2, cellInterpolationBasis2));

					assert(abs(weight1 + weight2 - 1) < Utils::NumericalZero);
				}
			});

		SparseMatrix Pi(finePb->HHO->nTotalFaceUnknowns, coarsePb->HHO->nElements * nCellUnknowns);
		parallelLoop.Fill(Pi);

		return Pi;
	}

	//-------------------------------------------------------------------------//
	//  Find polynomial faces which reconstruct the k+1 polynomial on the cell //
	//  and minimize the L2-norm                                               //
	//                     Min  1/2||x||^2_{L^2} = 1/2<x,M_F*x>                //
	//                     C*x = I                                             //
	//-------------------------------------------------------------------------//

	SparseMatrix GetGlobalMatrixFindFacesWhichReconstructCells(Diffusion_HHO<Dim>* problem)
	{
		if (!_useHigherOrderReconstruction)
			Utils::FatalError("This algo is only available if we reconstruct in degree k+1 on the cells.");

		// Re-init the faces because the mass matrix has been deleted after assembly for memory optimization,
		// but it is used in FindFacesPolyWhichReconstructOnTheCell(element)
		//problem->InitHHO_Faces();
		problem->InitHHO();

		int nFaceUnknowns = problem->HHO->nFaceUnknowns;
		int nCellUnknowns = problem->HHO->nReconstructUnknowns;

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, problem, nCellUnknowns, nFaceUnknowns](Element<Dim>* element, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* hhoElement = problem->HHOElement(element);
				DenseMatrix findFacesMatrix = FindFacesPolyWhichReconstructOnTheCell(hhoElement);

				for (auto face : element->Faces)
				{
					if (face->HasDirichletBC())
						continue;

					BigNumber elemGlobalNumber = element->Number;
					BigNumber faceGlobalNumber = face->Number;
					BigNumber faceLocalNumber = element->LocalNumberOf(face);

					double weight = Weight(element, face);
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, elemGlobalNumber*nCellUnknowns, weight*findFacesMatrix.block(faceLocalNumber*nFaceUnknowns, 0, nFaceUnknowns, nCellUnknowns));
				}
			});

		SparseMatrix M(problem->HHO->nTotalFaceUnknowns, problem->HHO->nElements * nCellUnknowns);
		parallelLoop.Fill(M);

		return M;
	}

	DenseMatrix FindFacesPolyWhichReconstructOnTheCell(Diff_HHOElement<Dim>* element)
	{
		//   | M_F  C^T | |x     |   |0|
		//   | C    0   | |lambda| = |I|

		auto HHO = element->HHO;

		// Assembly of the Lagrangian matrix
		int nBoundaryUnknowns = element->Faces.size() * HHO->nFaceUnknowns;
		DenseMatrix boundaryMassMatrix = DenseMatrix::Zero(nBoundaryUnknowns, nBoundaryUnknowns);
		for (auto face : element->Faces)
		{
			int localNumber = element->LocalNumberOf(face);
			boundaryMassMatrix.block(localNumber, localNumber, HHO->nFaceUnknowns, HHO->nFaceUnknowns) = face->MassMatrix();
		}

		DenseMatrix C = element->ReconstructionFromFacesMatrix();

		DenseMatrix lagrangianMatrix(nBoundaryUnknowns + C.rows(), nBoundaryUnknowns + C.rows());
		lagrangianMatrix.topLeftCorner(nBoundaryUnknowns, nBoundaryUnknowns) = boundaryMassMatrix;
		lagrangianMatrix.topRightCorner(C.cols(), C.rows()) = C.transpose();
		lagrangianMatrix.bottomLeftCorner(C.rows(), C.cols()) = C;
		lagrangianMatrix.bottomRightCorner(C.rows(), C.rows()) = DenseMatrix::Zero(C.rows(), C.rows());

		// Assembly of the right-hand side
		DenseMatrix rhs(HHO->nReconstructUnknowns + C.cols(), HHO->nReconstructUnknowns);
		rhs.topRows(C.cols()) = DenseMatrix::Zero(C.cols(), HHO->nReconstructUnknowns);
		rhs.bottomRows(HHO->nReconstructUnknowns) = DenseMatrix::Identity(HHO->nReconstructUnknowns, HHO->nReconstructUnknowns);

		// Solving
		Eigen::ColPivHouseholderQR<DenseMatrix> solver = lagrangianMatrix.colPivHouseholderQr();
		DenseMatrix solutionWithLagrangeCoeffs = solver.solve(rhs);
		DenseMatrix solution = solutionWithLagrangeCoeffs.topRows(nBoundaryUnknowns);

		return solution;
	}

	//---------------------//
	// Canonical injection //
	//---------------------//

	SparseMatrix GetGlobalCanonicalInjectionMatrixCoarseToFineElements()
	{
		Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		Mesh<Dim>* coarseMesh = coarsePb->_mesh;

		if (finePb->HHO->FaceBasis->GetDegree() < coarsePb->HHO->FaceBasis->GetDegree())
			Utils::FatalError("Canonical injection from coarse to fine elements only available if fine degree >= coarse degree.");

		int nCoarseUnknowns = _useHigherOrderReconstruction ? coarsePb->HHO->nReconstructUnknowns : coarsePb->HHO->nCellUnknowns;
		int nFineUnknowns = _useHigherOrderReconstruction ? finePb->HHO->nReconstructUnknowns : finePb->HHO->nCellUnknowns;

		ElementParallelLoop<Dim> parallelLoop(coarseMesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCoarseUnknowns * 4 * nFineUnknowns);

		parallelLoop.Execute([this, coarsePb, finePb, nCoarseUnknowns, nFineUnknowns](Element<Dim>* coarseElement, ParallelChunk<CoeffsChunk>* chunk)
			{
				DenseMatrix local_J_f_c = ComputeCanonicalInjectionMatrixCoarseToFine(coarseElement, coarsePb, finePb);
				for (auto fineElement : coarseElement->FinerElements)
				{
					BigNumber coarseElemGlobalNumber = coarseElement->Number;
					BigNumber fineElemGlobalNumber = fineElement->Number;
					BigNumber fineElemLocalNumber = coarseElement->LocalNumberOf(fineElement);

					chunk->Results.Coeffs.Add(fineElemGlobalNumber*nFineUnknowns, coarseElemGlobalNumber*nCoarseUnknowns, local_J_f_c.block(fineElemLocalNumber*nFineUnknowns, 0, nFineUnknowns, nCoarseUnknowns));
				}
			});

		SparseMatrix J_f_c(finePb->HHO->nElements * nFineUnknowns, coarsePb->HHO->nElements * nCoarseUnknowns);
		parallelLoop.Fill(J_f_c);
		return J_f_c;
	}

	DenseMatrix ComputeCanonicalInjectionMatrixCoarseToFine(Element<Dim>* coarseElement, Diffusion_HHO<Dim>* coarsePb, Diffusion_HHO<Dim>* finePb)
	{
		Diff_HHOElement<Dim>* hhoCoarseElem = coarsePb->HHOElement(coarseElement);
		FunctionalBasis<Dim>* coarseBasis = _useHigherOrderReconstruction ? hhoCoarseElem->ReconstructionBasis : hhoCoarseElem->CellBasis;
		int nCoarseUnknowns = coarseBasis->Size();
		int nFineUnknowns = _useHigherOrderReconstruction ? finePb->HHO->nReconstructUnknowns : finePb->HHO->nCellUnknowns;

		assert(nCoarseUnknowns <= nFineUnknowns);

		DenseMatrix J(nFineUnknowns * coarseElement->FinerElements.size(), nCoarseUnknowns);

		for (auto fineElement : coarseElement->FinerElements)
		{
			Diff_HHOElement<Dim>* hhoFineElem = finePb->HHOElement(fineElement);
			FunctionalBasis<Dim>* fineBasis = _useHigherOrderReconstruction ? hhoFineElem->ReconstructionBasis : hhoFineElem->CellBasis;

			DenseMatrix fineCoarseMass(fineBasis->Size(), coarseBasis->Size());
			for (BasisFunction<Dim>* finePhi : fineBasis->LocalFunctions)
			{
				for (BasisFunction<Dim>* coarsePhi : coarseBasis->LocalFunctions)
				{
					RefFunction functionToIntegrate = [coarseElement, fineElement, finePhi, coarsePhi](const RefPoint& fineRefPoint) {
						DomPoint domPoint = fineElement->ConvertToDomain(fineRefPoint);
						RefPoint coarseRefPoint = coarseElement->ConvertToReference(domPoint);
						return finePhi->Eval(fineRefPoint)*coarsePhi->Eval(coarseRefPoint);
					};

					int polynomialDegree = finePhi->GetDegree() + coarsePhi->GetDegree();
					double integral = fineElement->Integral(functionToIntegrate, polynomialDegree);
					fineCoarseMass(finePhi->LocalNumber, coarsePhi->LocalNumber) = integral;
				}
			}

			J.block(coarseElement->LocalNumberOf(fineElement)*nFineUnknowns, 0, nFineUnknowns, nCoarseUnknowns) = _useHigherOrderReconstruction ? hhoFineElem->SolveReconstructMassMatrix(fineCoarseMass) : hhoFineElem->SolveCellMassMatrix(fineCoarseMass);
		}

		return J;
	}


	//---------------------//
	//    L2 projection    //
	//---------------------//

	SparseMatrix GetGlobalL2ProjectionMatrixCoarseToFineElements()
	{
		Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		Mesh<Dim>* coarseMesh = coarsePb->_mesh;

		if (finePb->HHO->FaceBasis->GetDegree() < coarsePb->HHO->FaceBasis->GetDegree())
			Utils::FatalError("L2-proj from coarse to fine elements only available if fine degree > coarse degree.");

		int nCoarseUnknowns = _useHigherOrderReconstruction ? coarsePb->HHO->nReconstructUnknowns : coarsePb->HHO->nCellUnknowns;
		int nFineUnknowns = _useHigherOrderReconstruction ? finePb->HHO->nReconstructUnknowns : finePb->HHO->nCellUnknowns;

		ElementParallelLoop<Dim> parallelLoop(coarseMesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCoarseUnknowns * 6 * nFineUnknowns);

		parallelLoop.Execute([this, coarsePb, finePb, nCoarseUnknowns, nFineUnknowns](Element<Dim>* coarseElement, ParallelChunk<CoeffsChunk>* chunk)
			{
				DenseMatrix localL2Proj = ComputeL2ProjectionMatrixCoarseToFine(coarseElement, coarsePb, finePb);
				for (auto it = coarseElement->OverlappingFineElements.begin(); it != coarseElement->OverlappingFineElements.end(); it++)
				{
					Element<Dim>* fineElement = it->first;
					BigNumber coarseElemGlobalNumber = coarseElement->Number;
					BigNumber fineElemGlobalNumber = fineElement->Number;
					BigNumber fineElemLocalNumber = coarseElement->LocalNumberOfOverlapping(fineElement);

					chunk->Results.Coeffs.Add(fineElemGlobalNumber*nFineUnknowns, coarseElemGlobalNumber*nCoarseUnknowns, localL2Proj.block(fineElemLocalNumber*nFineUnknowns, 0, nFineUnknowns, nCoarseUnknowns));
				}
			});

		SparseMatrix L2Proj(finePb->HHO->nElements * nFineUnknowns, coarsePb->HHO->nElements * nCoarseUnknowns);
		parallelLoop.Fill(L2Proj);
		return L2Proj;
	}



	DenseMatrix ComputeL2ProjectionMatrixCoarseToFine(Element<Dim>* coarseElement, Diffusion_HHO<Dim>* coarsePb, Diffusion_HHO<Dim>* finePb)
	{
		assert(!coarseElement->OverlappingFineElements.empty());

		Diff_HHOElement<Dim>* hhoCoarseElem = coarsePb->HHOElement(coarseElement);
		FunctionalBasis<Dim>* coarseBasis = _useHigherOrderReconstruction ? hhoCoarseElem->ReconstructionBasis : hhoCoarseElem->CellBasis;
		int nCoarseUnknowns = coarseBasis->Size();
		int nFineUnknowns = _useHigherOrderReconstruction ? finePb->HHO->nReconstructUnknowns : finePb->HHO->nCellUnknowns;

		DenseMatrix L2Proj(nFineUnknowns * coarseElement->OverlappingFineElements.size(), nCoarseUnknowns);

		for (auto it = coarseElement->OverlappingFineElements.begin(); it != coarseElement->OverlappingFineElements.end(); it++)
		{
			Element<Dim>* fineElement = it->first;
			if (!fineElement->IsInSamePhysicalPartAs(coarseElement))
			{
				MatlabScript s;
				s.OpenFigure();
				string phypart = "no physical part";
				if (coarseElement->PhysicalPart)
					phypart = coarseElement->PhysicalPart->Name;
				s.Comment("---------------- Coarse element overlapped by the fine (phy part = " + phypart + ")");
				coarseElement->ExportToMatlab("b");
				phypart = "no physical part";
				if (coarseElement->PhysicalPart)
					phypart = fineElement->CoarserElement->PhysicalPart->Name;
				s.Comment("---------------- Coarse element associated to the fine (phy part = " + phypart + ")");
				fineElement->CoarserElement->ExportToMatlab("y");
				phypart = "no physical part";
				if (coarseElement->PhysicalPart)
					phypart = fineElement->PhysicalPart->Name;
				s.Comment("---------------- Fine element (phy part = " + phypart + ")");
				fineElement->ExportToMatlab("r");
				assert(false);
				Utils::FatalError("This coarse element is declared overlapped by a fine one that is not in the same physical part. The coarsening/refinement strategy must prevent that.");
			}

			Diff_HHOElement<Dim>* hhoFineElem = finePb->HHOElement(fineElement);
			FunctionalBasis<Dim>* fineBasis = _useHigherOrderReconstruction ? hhoFineElem->ReconstructionBasis : hhoFineElem->CellBasis;

			DenseMatrix fineCoarseMass(nFineUnknowns, nCoarseUnknowns);

			vector<PhysicalShape<Dim>*> intersectionCoarseFine = it->second;

			for (BasisFunction<Dim>* finePhi : fineBasis->LocalFunctions)
			{
				for (BasisFunction<Dim>* coarsePhi : coarseBasis->LocalFunctions)
				{
					double integral = 0;
					int degree = finePhi->GetDegree() + coarsePhi->GetDegree();
					for (PhysicalShape<Dim>* intersection : intersectionCoarseFine)
					{
						RefFunction finePhiCoarsePhi = [coarseElement, fineElement, intersection, finePhi, coarsePhi](const RefPoint& intersectionRefPoint) {
							DomPoint domPoint = intersection->ConvertToDomainAndSaveResult(intersectionRefPoint, true);
							RefPoint fineRefPoint = fineElement->ConvertToReference(domPoint);
							RefPoint coarseRefPoint = coarseElement->ConvertToReference(domPoint);
							return finePhi->Eval(fineRefPoint)*coarsePhi->Eval(coarseRefPoint);
						};
						integral += intersection->Integral(finePhiCoarsePhi, degree);
					}
					fineCoarseMass(finePhi->LocalNumber, coarsePhi->LocalNumber) = integral;
				}
			}

			L2Proj.block(coarseElement->LocalNumberOfOverlapping(fineElement)*nFineUnknowns, 0, nFineUnknowns, nCoarseUnknowns) = _useHigherOrderReconstruction ? hhoFineElem->SolveReconstructMassMatrix(fineCoarseMass) : hhoFineElem->SolveCellMassMatrix(fineCoarseMass);
		}

		return L2Proj;
	}

	//-----------------------------//
	// Canonical injection (faces) //
	//-----------------------------//

	SparseMatrix GetGlobalCanonicalInjectionMatrixCoarseToFineFaces()
	{
		Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		Mesh<Dim>* coarseMesh = coarsePb->_mesh;

		if (finePb->HHO->FaceBasis->GetDegree() != coarsePb->HHO->FaceBasis->GetDegree())
			Utils::FatalError("GetGlobalCanonicalInjectionMatrixCoarseToFineFaces() is not implemented when coarse and fine bases don't have the same degree.");

		int nFaceUnknowns = finePb->HHO->FaceBasis->Size();

		FaceParallelLoop<Dim> parallelLoop(coarseMesh->Faces);
		parallelLoop.ReserveChunkCoeffsSize(nFaceUnknowns * 2 * nFaceUnknowns);

		parallelLoop.Execute([this, nFaceUnknowns, coarsePb, finePb](Face<Dim>* coarseFace, ParallelChunk<CoeffsChunk>* chunk)
			{
				if (coarseFace->HasDirichletBC())
					return;

				DenseMatrix local_J_f_c = this->ComputeCanonicalInjectionMatrixCoarseToFine(coarseFace, coarsePb, finePb);
				for (auto fineFace : coarseFace->FinerFaces)
				{
					BigNumber coarseFaceGlobalNumber = coarseFace->Number;
					BigNumber fineFaceGlobalNumber = fineFace->Number;
					BigNumber fineFaceLocalNumber = coarseFace->LocalNumberOf(fineFace);

					chunk->Results.Coeffs.Add(fineFaceGlobalNumber*nFaceUnknowns, coarseFaceGlobalNumber*nFaceUnknowns, local_J_f_c.block(fineFaceLocalNumber*nFaceUnknowns, 0, nFaceUnknowns, nFaceUnknowns));
				}
			});

		SparseMatrix J_f_c(finePb->HHO->nInteriorAndNeumannFaces * nFaceUnknowns, coarsePb->HHO->nInteriorAndNeumannFaces * nFaceUnknowns);
		parallelLoop.Fill(J_f_c);
		return J_f_c;
	}



	DenseMatrix ComputeCanonicalInjectionMatrixCoarseToFine(Face<Dim>* coarseFace, Diffusion_HHO<Dim>* coarseProblem, Diffusion_HHO<Dim>* fineProblem)
	{
		int coarseFaceUnknowns = coarseProblem->HHO->nFaceUnknowns;
		int fineFaceUnknowns = fineProblem->HHO->nFaceUnknowns;
		assert(coarseFaceUnknowns <= fineFaceUnknowns);

		DenseMatrix J(fineFaceUnknowns * coarseFace->FinerFaces.size(), coarseFaceUnknowns);

		Diff_HHOFace<Dim>* coarseHHOFace = coarseProblem->HHOFace(coarseFace);

		for (auto fineFace : coarseFace->FinerFaces)
		{
			Diff_HHOFace<Dim>* fineHHOFace = fineProblem->HHOFace(fineFace);

			DenseMatrix fineCoarseMass(fineFaceUnknowns, coarseFaceUnknowns);
			for (BasisFunction<Dim - 1>* finePhi : fineHHOFace->Basis->LocalFunctions)
			{
				for (BasisFunction<Dim - 1>* coarsePhi : coarseHHOFace->Basis->LocalFunctions)
				{
					RefFunction functionToIntegrate = [coarseFace, fineFace, finePhi, coarsePhi](const RefPoint& fineRefPoint) {
						DomPoint domPoint = fineFace->ConvertToDomain(fineRefPoint);
						RefPoint coarseRefPoint = coarseFace->ConvertToReference(domPoint);
						return finePhi->Eval(fineRefPoint)*coarsePhi->Eval(coarseRefPoint);
					};

					int polynomialDegree = finePhi->GetDegree() + coarsePhi->GetDegree();
					double integral = fineFace->Integral(functionToIntegrate, polynomialDegree);
					fineCoarseMass(finePhi->LocalNumber, coarsePhi->LocalNumber) = integral;
				}
			}
			J.block(coarseFace->LocalNumberOf(fineFace)*fineFaceUnknowns, 0, fineFaceUnknowns, coarseFaceUnknowns) = fineHHOFace->SolveMassMatrix(fineCoarseMass);
		}

		return J;
	}

	//-------------------------------------------//
	//  Used by the algorithm from Wildey et al. //
	//-------------------------------------------//

	DenseMatrix StaticallyCondenseInteriorFinerFaces(Element<Dim>* coarseElem, const SparseMatrix& fineA, Diffusion_HHO<Dim>* coarsePb, Diffusion_HHO<Dim>* finePb)
	{
		int nFaceUnknowns = coarsePb->HHO->nFaceUnknowns;

		int nFineInteriorFaces = coarseElem->FinerFacesRemoved.size();
		int nFineBoundaryFaces = 0;
		for (Face<Dim>* cf : coarseElem->Faces)
			nFineBoundaryFaces += cf->FinerFaces.size();

		DenseMatrix Aii = DenseMatrix::Zero(nFineInteriorFaces*nFaceUnknowns, nFineInteriorFaces*nFaceUnknowns);
		DenseMatrix Aib = DenseMatrix::Zero(nFineInteriorFaces*nFaceUnknowns, nFineBoundaryFaces*nFaceUnknowns);

		for (int i = 0; i < coarseElem->FinerFacesRemoved.size(); i++)
		{
			// Aii
			Face<Dim>* fi = coarseElem->FinerFacesRemoved[i];
			for (int j = 0; j < coarseElem->FinerFacesRemoved.size(); j++)
			{
				Face<Dim>* fj = coarseElem->FinerFacesRemoved[j];
				Aii.block(i*nFaceUnknowns, j*nFaceUnknowns, nFaceUnknowns, nFaceUnknowns) = fineA.block(fi->Number*nFaceUnknowns, fj->Number*nFaceUnknowns, nFaceUnknowns, nFaceUnknowns);
			}
			// Aib
			BigNumber fineFaceLocalNumberInCoarseElem = 0;
			for (Face<Dim>* cf : coarseElem->Faces)
			{
				if (cf->IsDomainBoundary)
					continue;

				for (int j = 0; j < cf->FinerFaces.size(); j++)
				{
					Face<Dim>* fj = cf->FinerFaces[j];
					Aib.block(i*nFaceUnknowns, fineFaceLocalNumberInCoarseElem*nFaceUnknowns, nFaceUnknowns, nFaceUnknowns) = fineA.block(fi->Number*nFaceUnknowns, fj->Number*nFaceUnknowns, nFaceUnknowns, nFaceUnknowns);
					fineFaceLocalNumberInCoarseElem++;
				}
			}
		}


		int nCoarseFaces = 0;
		for (Face<Dim>* cf : coarseElem->Faces)
			nCoarseFaces++;

		DenseMatrix J = DenseMatrix::Zero(nFineBoundaryFaces*nFaceUnknowns, nCoarseFaces*nFaceUnknowns);
		BigNumber fineFaceLocalNumberInCoarseElem = 0;
		for (auto coarseFace : coarseElem->Faces)
		{
			if (coarseFace->IsDomainBoundary)
				continue;

			DenseMatrix local_J_f_c = ComputeCanonicalInjectionMatrixCoarseToFine(coarseFace, coarsePb, finePb);
			for (auto fineFace : coarseFace->FinerFaces)
			{
				BigNumber fineFaceLocalNumberInCoarseFace = coarseFace->LocalNumberOf(fineFace);
				BigNumber coarseFaceLocalNumberInCoarseElem = coarseElem->LocalNumberOf(coarseFace);

				J.block(fineFaceLocalNumberInCoarseElem*nFaceUnknowns, coarseFaceLocalNumberInCoarseElem*nFaceUnknowns, nFaceUnknowns, nFaceUnknowns) = local_J_f_c.block(fineFaceLocalNumberInCoarseFace*nFaceUnknowns, 0, nFaceUnknowns, nFaceUnknowns);
				fineFaceLocalNumberInCoarseElem++;
			}
		}

		DenseMatrix resolveCondensedFinerFacesFromFineBoundary = -Aii.llt().solve(Aib);
		DenseMatrix resolveCondensedFinerFacesFromCoarseBoundary = resolveCondensedFinerFacesFromFineBoundary * J;

		return resolveCondensedFinerFacesFromCoarseBoundary;
	}


public:
	~LevelForHHO()
	{
		if (!this->IsFinestLevel())
			delete _problem;
	}
};