#pragma once
#include "Multigrid.h"
#include "../HHO/Diffusion_HHO.h"
#include "../Utils/ElementParallelLoop.h"
using namespace std;

template <int Dim>
class LevelForHHO : public Level
{
private:
	Prolongation _prolongationCode = Prolongation::CellInterp_Trace;
	FunctionalBasis<Dim>* _cellInterpolationBasis;
	string _weightCode = "k";
public:
	Diffusion_HHO<Dim>* _problem;

	LevelForHHO(int number, Diffusion_HHO<Dim>* problem, Prolongation prolongationCode, FunctionalBasis<Dim>* cellInterpolationBasis, string weightCode)
		: Level(number)
	{
		this->_problem = problem;
		this->_prolongationCode = prolongationCode;
		this->_cellInterpolationBasis = cellInterpolationBasis;
		this->_weightCode = weightCode;
	}

	BigNumber NUnknowns() override
	{
		return _problem->HHO->nTotalFaceUnknowns;
	}

	void CoarsenMesh(CoarseningStrategy coarseningStgy, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) override
	{
		if (Utils::IsRefinementStrategy(coarseningStgy) && _problem->_mesh->CoarseMesh == nullptr)
		{
			noCoarserMeshProvided = true;
			return;
		}
		cout << "\tCoarsening mesh" << endl;
		_problem->_mesh->CoarsenMesh(coarseningStgy);
		if (_problem->_mesh->CoarseMesh->InteriorFaces.size() == 0)
			coarsestPossibleMeshReached = true;
	}

	void ExportVector(const Vector& v, string suffix, int levelNumber) override
	{
		if (!this->IsFinestLevel())
			this->FinerLevel->ExportVector(v, suffix, levelNumber);
		else
			this->_problem->ExportVector(v, "level" + to_string(levelNumber) + "_" + suffix);
	}

	void ExportMatrix(const SparseMatrix& M, string suffix, int levelNumber) override
	{
		if (!this->IsFinestLevel())
			this->FinerLevel->ExportMatrix(M, suffix, levelNumber);
		else
			this->_problem->ExportMatrix(M, "level" + to_string(levelNumber) + "_" + suffix);
	}

	void ExportMeshToMatlab(Mesh<Dim>* levelMesh, int levelNumber)
	{
		if (!this->IsFinestLevel())
			dynamic_cast<LevelForHHO<Dim>*>(this->FinerLevel)->ExportMeshToMatlab(levelMesh, levelNumber);
		else
		{
			string filePath = this->_problem->GetFilePath("level" + to_string(levelNumber) + "_faces");
			levelMesh->ExportFacesToMatlab(filePath);
		}
	}

private:
	void SetupDiscretizedOperator() override 
	{
		this->OperatorMatrix = &this->_problem->A;
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
		if (_prolongationCode == Prolongation::FaceInject)
			cout << "\t\tMesh                : " << this->_problem->_mesh->Faces.size() << " faces" << endl;
		else
		{
			cout << "\t\tMesh                : " << this->_problem->_mesh->Elements.size() << " elements, regularity = ";
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
			this->ExportMeshToMatlab(this->_problem->_mesh, this->Number);

		if (!IsCoarsestLevel())
		{
			Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

			if (this->UseGalerkinOperator)
				coarsePb->InitHHO();
			else
			{
				ActionsArguments actions;
				actions.LogAssembly = false;
				actions.AssembleRightHandSide = false;
				coarsePb->Assemble(actions);
			}
		}
	}

	void SetupProlongation() override
	{
		Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

		if (_prolongationCode == Prolongation::CellInterp_Trace)
		{
			//----------------------------------------------------------//
			//                       Default method                     //
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
		else if (_prolongationCode == Prolongation::CellInterp_Inject_Trace)
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
		else if (_prolongationCode == Prolongation::CellInterp_L2proj_Trace)
		{
			//---------------------------------------------------------------------//
			//                      Non-nested variant                             //
			// Step 1: Interpolation from coarse faces to coarse cells.            //
			// Step 2: Instead of the canonical injection which is now impossible, //
			//         L2-projection onto the fine cells.                          //
			// Step 3: Trace on the fine faces.                                    //
			//---------------------------------------------------------------------//

			CoarseningStrategy stgy = coarsePb->_mesh->ComesFrom.CS;
			if (stgy == CoarseningStrategy::None)
				stgy = finePb->_mesh->ComesFrom.CS;

			// CheckIfFullyEmbeddedInCoarseElement
			ElementParallelLoop<Dim> parallelLoop(finePb->_mesh->Elements);
			parallelLoop.Execute([stgy](Element<Dim>* fe, ParallelChunk<CoeffsChunk>* chunk)
				{
					fe->CheckIfFullyEmbeddedInCoarseElement(stgy);
				});

			// SetOverlappingFineElements
			ElementParallelLoop<Dim> parallelLoop2(coarsePb->_mesh->Elements);
			parallelLoop2.Execute([stgy](Element<Dim>* ce, ParallelChunk<CoeffsChunk>* chunk)
				{
					ce->SetOverlappingFineElements(stgy);
					ce->InitOverlappingElementsLocalNumbering();
				});

			SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
			SparseMatrix J_f_c = GetGlobalL2ProjectionMatrixCoarseToFineElements();
			SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCellsToFaces(finePb);

			if (ExportComponents)
			{
				Level::ExportMatrix(I_c, "I_c");
				Level::ExportMatrix(J_f_c, "J_f_c");
				Level::ExportMatrix(Pi_f, "Pi_f");
			}

			P = Pi_f * J_f_c * I_c;
		}
		else if (_prolongationCode == Prolongation::CellInterp_InjectAndTrace)
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

			parallelLoop.Execute([this, P_algo1, J_faces, nFaceUnknowns](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
				{
					if (f->HasDirichletBC())
						return;

					Diff_HHOFace<Dim>* fineFace = dynamic_cast<Diff_HHOFace<Dim>*>(f);
					if (f->IsRemovedOnCoarserGrid)
						chunk->Results.Coeffs.CopyRows(fineFace->Number*nFaceUnknowns, nFaceUnknowns, P_algo1);
					else
						chunk->Results.Coeffs.CopyRows(fineFace->Number*nFaceUnknowns, nFaceUnknowns, J_faces);
				});

			P = SparseMatrix(finePb->HHO->nInteriorAndNeumannFaces * nFaceUnknowns, coarsePb->HHO->nInteriorAndNeumannFaces * nFaceUnknowns);
			parallelLoop.Fill(P);
		}
		else if (_prolongationCode == Prolongation::CellInterp_Inject_Adjoint)
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
		else if (_prolongationCode == Prolongation::Wildey)
		{
			//-------------------------------------------------------------------------------------------------//
			//                                           Wildey et al.                                         //
			// The coarse level is built by static condensation of the fine faces interior to coarse elements. //
			// The prolongation solves those condensed unknowns.                                               //
			//-------------------------------------------------------------------------------------------------//

			int nFaceUnknowns = finePb->HHO->nFaceUnknowns;

			ElementParallelLoop<Dim> parallelLoop(coarsePb->_mesh->Elements);

			parallelLoop.Execute([this, finePb, nFaceUnknowns](Element<Dim>* ce, ParallelChunk<CoeffsChunk>* chunk)
				{
					Diff_HHOElement<Dim>* coarseElem = dynamic_cast<Diff_HHOElement<Dim>*>(ce);
					DenseMatrix resolveCondensedFinerFacesFromCoarseBoundary = coarseElem->StaticallyCondenseInteriorFinerFaces(*this->OperatorMatrix);

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

					for (Face<Dim>* cf : coarseElem->Faces)
					{
						if (cf->HasDirichletBC())
							continue;

						Diff_HHOFace<Dim>* coarseFace = dynamic_cast<Diff_HHOFace<Dim>*>(cf);
						DenseMatrix local_J_f_c = coarseFace->ComputeCanonicalInjectionMatrixCoarseToFine(finePb->HHO->FaceBasis);
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
		else if (_prolongationCode == Prolongation::FaceInject)
		{
			//----------------------------------------------------------------//
			// Canonical injection from coarse faces to fine faces.           //
			// Implemented to be used with the face coarsening (option -cs f) //
			//----------------------------------------------------------------//

			SparseMatrix J_faces = GetGlobalCanonicalInjectionMatrixCoarseToFineFaces();
			P = J_faces;
		}
		else
			Utils::FatalError("Unknown prolongationCode.");
	}

	void SetupRestriction() override
	{
		R = RestrictionScalingFactor() * P.transpose();
	}

	void OnEndSetup() override
	{
		// If finest level, delete everything you don't need to reconstruct the solution at the end
		if (IsFinestLevel())
		{
			ElementParallelLoop<Dim> parallelLoopE(_problem->_mesh->Elements);
			parallelLoopE.Execute([](Element<Dim>* element)
				{
					Diff_HHOElement<Dim>* e = dynamic_cast<Diff_HHOElement<Dim>*>(element);
					e->DeleteUselessMatricesAfterMultigridSetup();
				});

			FaceParallelLoop<Dim> parallelLoopF(_problem->_mesh->Faces);
			parallelLoopF.Execute([](Face<Dim>* face)
				{
					Diff_HHOFace<Dim>* f = dynamic_cast<Diff_HHOFace<Dim>*>(face);
					f->DeleteUselessMatricesAfterMultigridSetup();
				});
		}

		// On the coarse levels, delete the whole mesh
		_problem->_mesh->CoarseMesh = nullptr;
		if (!IsFinestLevel())
			delete _problem->_mesh;
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
		DomFunction constantOne = [](DomPoint p) { return 1; };
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
		DomFunction constantOne = [](DomPoint p) { return 1; };
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
		else if (_weightCode.compare("k") == 0)
		{
			auto n = element->OuterNormalVector(face);
			double k = (element->DiffTensor() * n).dot(n);

			auto n1 = face->Element1->OuterNormalVector(face);
			double k1 = (face->Element1->DiffTensor() * n1).dot(n1);

			auto n2 = face->Element2->OuterNormalVector(face);
			double k2 = (face->Element2->DiffTensor() * n2).dot(n2);

			return k / (k1 + k2);
		}
		else if (_weightCode.compare("a") == 0)
			return 0.5;
		else
			Utils::FatalError("Unknown weight code.");
		assert(false);
	}

	SparseMatrix GetGlobalInterpolationMatrixFromFacesToCells(Diffusion_HHO<Dim>* problem)
	{
		int nCellUnknowns = _cellInterpolationBasis->Size();
		int nFaceUnknowns = problem->HHO->nFaceUnknowns;

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nFaceUnknowns, nCellUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* element = dynamic_cast<Diff_HHOElement<Dim>*>(e);

				DenseMatrix cellInterpMatrix;
				if (_cellInterpolationBasis == _problem->HHO->ReconstructionBasis)
					cellInterpMatrix = element->ReconstructionFromFacesMatrix();
				else if (_cellInterpolationBasis == _problem->HHO->CellBasis)
					cellInterpMatrix = element->SolveCellUnknownsMatrix();
				else
					assert(false);

				for (auto f : element->Faces)
				{
					if (f->HasDirichletBC())
						continue;

					Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);

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

		parallelLoop.Execute([this, nFaceUnknowns, nCellUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* element = dynamic_cast<Diff_HHOElement<Dim>*>(e);
				DenseMatrix reconstructMatrix = element->SolveCellUnknownsMatrix();

				for (auto f : element->Faces)
				{
					if (f->HasDirichletBC())
						continue;

					Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);

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
		int nFaceUnknowns = problem->HHO->nFaceUnknowns;
		int nCellUnknowns = _cellInterpolationBasis->Size();

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nCellUnknowns, nFaceUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* element = dynamic_cast<Diff_HHOElement<Dim>*>(e);

				for (auto f : element->Faces)
				{
					if (f->HasDirichletBC())
						continue;

					Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);

					BigNumber elemGlobalNumber = element->Number;
					BigNumber faceGlobalNumber = face->Number;

					double weight = Weight(element, face);
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, elemGlobalNumber*nCellUnknowns, weight*face->Trace(element, _cellInterpolationBasis));
				}
			});

		SparseMatrix Pi(problem->HHO->nTotalFaceUnknowns, problem->HHO->nElements * nCellUnknowns);
		parallelLoop.Fill(Pi);

		return Pi;
	}

	SparseMatrix GetGlobalProjectorMatrixFromCoarseCellsToFineFaces()
	{
		Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

		int nFaceUnknowns = finePb->HHO->nFaceUnknowns;
		int nCellUnknowns = _cellInterpolationBasis->Size();

		FaceParallelLoop<Dim> parallelLoop(finePb->_mesh->Faces);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nCellUnknowns, nFaceUnknowns](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				if (f->HasDirichletBC())
					return;

				Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);
				BigNumber faceGlobalNumber = face->Number;
				if (face->IsDomainBoundary || face->IsRemovedOnCoarserGrid)
				{
					Element<Dim>* coarseElem = face->Element1->CoarserElement;
					BigNumber coarseElemGlobalNumber = coarseElem->Number;

					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, coarseElemGlobalNumber*nCellUnknowns, face->Trace(coarseElem, _cellInterpolationBasis));
				}
				else
				{
					Element<Dim>* coarseElem1 = face->Element1->CoarserElement;
					Element<Dim>* coarseElem2 = face->Element2->CoarserElement;
					BigNumber coarseElem1GlobalNumber = coarseElem1->Number;
					BigNumber coarseElem2GlobalNumber = coarseElem2->Number;

					//double weight1 = Weight(coarseElem1, face->CoarseFace); // Careful with using face->CoarseFace when the mesh isn't nested...
					double weight1 = Weight(face->Element1, face);
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, coarseElem1GlobalNumber*nCellUnknowns, weight1*face->Trace(coarseElem1, _cellInterpolationBasis));

					//double weight2 = Weight(coarseElem2, face->CoarseFace);
					double weight2 = Weight(face->Element2, face);
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, coarseElem2GlobalNumber*nCellUnknowns, weight2*face->Trace(coarseElem2, _cellInterpolationBasis));

					assert(abs(weight1 + weight2 - 1) < 1e-12);
				}
			});

		SparseMatrix Pi(finePb->HHO->nTotalFaceUnknowns, coarsePb->HHO->nElements * nCellUnknowns);
		parallelLoop.Fill(Pi);

		return Pi;
	}

	// Min 1/2 ||  ||^2
	SparseMatrix GetGlobalMatrixFindFacesWhichReconstructCells(Diffusion_HHO<Dim>* problem)
	{
		if (_cellInterpolationBasis != _problem->HHO->ReconstructionBasis)
			assert(false && "This algo is only available if we reconstruct in degree k+1 on the cells.");

		int nFaceUnknowns = problem->HHO->nFaceUnknowns;
		int nCellUnknowns = _cellInterpolationBasis->Size();

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nCellUnknowns, nFaceUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* element = dynamic_cast<Diff_HHOElement<Dim>*>(e);

				DenseMatrix findFacesMatrix = element->FindFacesPolyWhichReconstructOnTheCell();

				for (auto f : element->Faces)
				{
					if (f->HasDirichletBC())
						continue;

					Diff_HHOFace<Dim>* face = dynamic_cast<Diff_HHOFace<Dim>*>(f);

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

	SparseMatrix GetGlobalCanonicalInjectionMatrixCoarseToFineElements()
	{
		Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		Mesh<Dim>* coarseMesh = coarsePb->_mesh;

		int nCellUnknowns = _cellInterpolationBasis->Size();

		ElementParallelLoop<Dim> parallelLoop(coarseMesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nCellUnknowns);

		parallelLoop.Execute([this, nCellUnknowns](Element<Dim>* ce, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* coarseElement = dynamic_cast<Diff_HHOElement<Dim>*>(ce);

				DenseMatrix local_J_f_c = coarseElement->ComputeCanonicalInjectionMatrixCoarseToFine(_cellInterpolationBasis);
				for (auto fineElement : coarseElement->FinerElements)
				{
					BigNumber coarseElemGlobalNumber = coarseElement->Number;
					BigNumber fineElemGlobalNumber = fineElement->Number;
					BigNumber fineElemLocalNumber = coarseElement->LocalNumberOf(fineElement);

					chunk->Results.Coeffs.Add(fineElemGlobalNumber*nCellUnknowns, coarseElemGlobalNumber*nCellUnknowns, local_J_f_c.block(fineElemLocalNumber*nCellUnknowns, 0, nCellUnknowns, nCellUnknowns));
				}
			});

		SparseMatrix J_f_c(finePb->HHO->nElements * nCellUnknowns, coarsePb->HHO->nElements * nCellUnknowns);
		parallelLoop.Fill(J_f_c);
		return J_f_c;
	}

	SparseMatrix GetGlobalL2ProjectionMatrixCoarseToFineElements()
	{
		Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		Mesh<Dim>* coarseMesh = coarsePb->_mesh;

		int nCellUnknowns = _cellInterpolationBasis->Size();

		ElementParallelLoop<Dim> parallelLoop(coarseMesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nCellUnknowns);

		parallelLoop.Execute([this, nCellUnknowns](Element<Dim>* ce, ParallelChunk<CoeffsChunk>* chunk)
			{
				Diff_HHOElement<Dim>* coarseElement = dynamic_cast<Diff_HHOElement<Dim>*>(ce);

				DenseMatrix local_J_f_c = coarseElement->ComputeL2ProjectionMatrixCoarseToFine(_cellInterpolationBasis);
				for (auto fineElement : coarseElement->OverlappingFineElements)
				{
					BigNumber coarseElemGlobalNumber = coarseElement->Number;
					BigNumber fineElemGlobalNumber = fineElement->Number;
					BigNumber fineElemLocalNumber = coarseElement->LocalNumberOfOverlapping(fineElement);

					chunk->Results.Coeffs.Add(fineElemGlobalNumber*nCellUnknowns, coarseElemGlobalNumber*nCellUnknowns, local_J_f_c.block(fineElemLocalNumber*nCellUnknowns, 0, nCellUnknowns, nCellUnknowns));
				}
			});

		SparseMatrix J_f_c(finePb->HHO->nElements * nCellUnknowns, coarsePb->HHO->nElements * nCellUnknowns);
		parallelLoop.Fill(J_f_c);
		return J_f_c;
	}

	SparseMatrix GetGlobalCanonicalInjectionMatrixCoarseToFineFaces()
	{
		Diffusion_HHO<Dim>* finePb = this->_problem;
		Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		Mesh<Dim>* coarseMesh = coarsePb->_mesh;

		FunctionalBasis<Dim - 1>* faceBasis = finePb->HHO->FaceBasis;
		int nFaceUnknowns = faceBasis->Size();

		FaceParallelLoop<Dim> parallelLoop(coarseMesh->Faces);
		parallelLoop.ReserveChunkCoeffsSize(nFaceUnknowns * 2 * nFaceUnknowns);

		parallelLoop.Execute([this, faceBasis, nFaceUnknowns](Face<Dim>* cf, ParallelChunk<CoeffsChunk>* chunk)
			{
				if (cf->HasDirichletBC())
					return;

				Diff_HHOFace<Dim>* coarseFace = dynamic_cast<Diff_HHOFace<Dim>*>(cf);

				DenseMatrix local_J_f_c = coarseFace->ComputeCanonicalInjectionMatrixCoarseToFine(faceBasis);
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

public:
	~LevelForHHO()
	{
		if (!this->IsFinestLevel())
			delete _problem;
	}
};

template <int Dim>
class MultigridForHHO : public Multigrid
{
private:
	Diffusion_HHO<Dim>* _problem;
	Prolongation _prolongationCode = Prolongation::CellInterp_Trace;
	FunctionalBasis<Dim>* _cellInterpolationBasis;
	string _weightCode = "k";
public:

	MultigridForHHO(Diffusion_HHO<Dim>* problem, Prolongation prolongationCode, FunctionalBasis<Dim>* cellInterpolationBasis, string weightCode)
		: MultigridForHHO(problem, prolongationCode, cellInterpolationBasis, weightCode, 0)
	{}

	MultigridForHHO(Diffusion_HHO<Dim>* problem, Prolongation prolongationCode, FunctionalBasis<Dim>* cellInterpolationBasis, string weightCode, int nLevels)
		: Multigrid(nLevels)
	{
		this->_problem = problem;
		this->_prolongationCode = prolongationCode;
		this->_cellInterpolationBasis = cellInterpolationBasis;
		this->_weightCode = weightCode;
		this->BlockSizeForBlockSmoothers = problem->HHO->nFaceUnknowns;
		this->_fineLevel = new LevelForHHO<Dim>(0, problem, _prolongationCode, _cellInterpolationBasis, _weightCode);
	}

	void BeginSerialize(ostream& os) const override
	{
		os << "MultigridForHHO" << endl;
		os << "\t" << "Prolongation       : ";
		if (_prolongationCode == Prolongation::CellInterp_Trace)
			os << "coarse cell interpolation + trace on fine faces ";
		else if (_prolongationCode == Prolongation::CellInterp_Inject_Trace)
			os << "coarse cell interpolation + injection coarse to fine cells + trace on fine faces ";
		else if (_prolongationCode == Prolongation::CellInterp_L2proj_Trace)
			os << "coarse cell interpolation + L2-projection to fine cells + trace on fine faces ";
		else if (_prolongationCode == Prolongation::CellInterp_InjectAndTrace)
			os << "injection for common faces, and coarse cell interpolation + trace for the other ";
		else if (_prolongationCode == Prolongation::CellInterp_Inject_Adjoint)
			os << "coarse cell interpolation + injection coarse to fine cells + adjoint of cell interpolation ";
		else if (_prolongationCode == Prolongation::Wildey)
			os << "Wildey et al. ";
		else if (_prolongationCode == Prolongation::FaceInject)
			os << "injection coarse to fine faces ";
		os << "[-prolong " << (unsigned)_prolongationCode << "]" << endl;

		if (_prolongationCode != Prolongation::FaceInject)
		{
			os << "\t" << "Cell interpolation : ";
			if (_cellInterpolationBasis == _problem->HHO->CellBasis)
				os << "k [-cell-interp 2]";
			else if (_cellInterpolationBasis == _problem->HHO->ReconstructionBasis)
				os << "k+1 [-cell-interp 1]";
			os << endl;
		}

		os << "\t" << "Weighting          : ";
		if (_weightCode.compare("k") == 0)
			os << "proportional to the diffusion coefficient [-weight k]";
		else if (_weightCode.compare("a") == 0)
			os << "simple average [-weight a]";
		os << endl;
	}

	void EndSerialize(ostream& os) const override
	{
		if ((this->CoarseningStgy == CoarseningStrategy::AgglomerationCoarseningByFaceNeighbours
			|| this->CoarseningStgy == CoarseningStrategy::AgglomerationCoarseningByClosestCenter
			|| this->CoarseningStgy == CoarseningStrategy::AgglomerationCoarseningByClosestFace
			|| this->CoarseningStgy == CoarseningStrategy::AgglomerationCoarseningByLargestInterface
			|| this->CoarseningStgy == CoarseningStrategy::AgglomerationCoarseningByVertexNeighbours)
			&& _prolongationCode != Prolongation::CellInterp_L2proj_Trace)
		{
			os << endl;
			Utils::Warning(os, "The selected coarsening strategy generates non-nested meshes, while the selected prolongation operator is made for nested meshes. Option '-prolong " + to_string((unsigned)Prolongation::CellInterp_L2proj_Trace) + "' recommended.");
		}
	}

	Vector Solve(const Vector& b, string initialGuessCode) override
	{
		Vector initialGuess;
		if (initialGuessCode.compare("smooth") == 0)
		{
			DomFunction coarseErrorFunction = [](DomPoint p)
			{
				double x = p.X;
				double y = p.Y;
				//return sin(4 * M_PI * x)*sin(4 * M_PI * y);
				//return x * (1 - x) * y*(1 - y);
				return sin(M_PI * x)*sin(M_PI * y);
			};

			Level* coarseLevel = this->_fineLevel;
			while (coarseLevel->CoarserLevel)
				coarseLevel = coarseLevel->CoarserLevel;

			Diffusion_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(coarseLevel)->_problem;

			Vector coarseError = coarsePb->ProjectOnFaceDiscreteSpace(coarseErrorFunction);

			Level* fineLevel = coarseLevel;
			Vector prolongatedError = coarseError;
			while (fineLevel->FinerLevel)
			{
				fineLevel = fineLevel->FinerLevel;
				prolongatedError = fineLevel->Prolong(prolongatedError);
			}

			initialGuess = -prolongatedError;
		}
		else
			return Multigrid::Solve(b, initialGuessCode);

		return Multigrid::Solve(b, initialGuess);
	}
	
protected:
	Level* CreateCoarseLevel(Level* fineLevel) override
	{
		LevelForHHO<Dim>* hhoFineLevel = dynamic_cast<LevelForHHO<Dim>*>(fineLevel);
		Diffusion_HHO<Dim>* coarseProblem = hhoFineLevel->_problem->GetProblemOnCoarserMesh();
		LevelForHHO<Dim>* coarseLevel = new LevelForHHO<Dim>(fineLevel->Number + 1, coarseProblem, _prolongationCode, _cellInterpolationBasis, _weightCode);
		return coarseLevel;
	}
};