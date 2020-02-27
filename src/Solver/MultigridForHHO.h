#pragma once
#include "Multigrid.h"
#include "../HHO/Poisson_HHO.h"
#include "../Utils/ElementParallelLoop.h"
using namespace std;

template <int Dim>
class LevelForHHO : public Level
{
private:
	int _prolongationCode = 1;
	FunctionalBasis<Dim>* _cellInterpolationBasis;
public:
	Poisson_HHO<Dim>* _problem;

	LevelForHHO(int number, Poisson_HHO<Dim>* problem, int prolongationCode, FunctionalBasis<Dim>* cellInterpolationBasis)
		: Level(number)
	{
		this->_problem = problem;
		this->_prolongationCode = prolongationCode;
		this->_cellInterpolationBasis = cellInterpolationBasis;

		if (ExportMatrices)
			problem->ExportFaces("L" + to_string(number));
	}

	BigNumber NUnknowns() override
	{
		return _problem->HHO->nTotalFaceUnknowns;
	}

	void CoarsenMesh(CoarseningStrategy coarseningStgy, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) override
	{
		if (coarseningStgy == CoarseningStrategy::SplittingRefinement && _problem->_mesh->CoarseMesh == nullptr)
		{
			noCoarserMeshProvided = true;
			return;
		}
		_problem->_mesh->CoarsenMesh(coarseningStgy);
		if (_problem->_mesh->CoarseMesh->InteriorFaces.size() == 0)
			coarsestPossibleMeshReached = true;
	}

	void ExportVector(Vector& v, string suffix) override
	{
		this->_problem->ExportVector(v, suffix);
	}

private:
	void SetupDiscretizedOperator() override 
	{
		this->OperatorMatrix = this->_problem->A;
	}

	void OnStartSetup() override
	{
		if (_prolongationCode == 5)
			cout << "\t\tMesh                : " << this->_problem->_mesh->Faces.size() << " faces" << endl;
		else
			cout << "\t\tMesh                : " << this->_problem->_mesh->Elements.size() << " elements, regularity = " << this->_problem->_mesh->Regularity() << endl;

		if (!IsCoarsestLevel())
		{
			Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

			if (this->UseGalerkinOperator)
				coarsePb->InitHHO();
			else
				coarsePb->Assemble(Action::None);
		}
	}

	void SetupProlongation() override
	{
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

		if (_prolongationCode == 1)
		{
			SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
			//SparseMatrix J_f_c = GetGlobalCanonicalInjectionMatrixCoarseToFineElements();
			//SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCellsToFaces(finePb);
			SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCoarseCellsToFineFaces();

			if (ExportMatrices)
			{
				finePb->ExportMatrix(I_c, "I_c");
				//finePb->ExportMatrix(J_f_c, "J_f_c");
				finePb->ExportMatrix(Pi_f, "Pi_f");

				SparseMatrix SolveCellUnknowns = GetSolveCellUnknownsMatrix(coarsePb);
				finePb->ExportMatrix(SolveCellUnknowns, "SolveCellUnknowns");
			}

			//P = Pi_f * J_f_c * I_c;
			P = Pi_f * I_c;
		}
		else if (_prolongationCode == 2)
		{
			SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
			//SparseMatrix J_f_c = GetGlobalCanonicalInjectionMatrixCoarseToFineElements();
			//SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCellsToFaces(finePb);
			SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCoarseCellsToFineFaces();

			SparseMatrix J_faces = GetGlobalCanonicalInjectionMatrixCoarseToFineFaces();

			if (ExportMatrices)
			{
				finePb->ExportMatrix(I_c, "I_c");
				//finePb->ExportMatrix(J_f_c, "J_f_c");
				finePb->ExportMatrix(Pi_f, "Pi_f");

				finePb->ExportMatrix(J_faces, "J_faces");
			}

			//SparseMatrix P_algo1 = Pi_f * J_f_c * I_c;
			SparseMatrix P_algo1 = Pi_f * I_c;
			
			int nFaceUnknowns = finePb->HHO->nFaceUnknowns;

			FaceParallelLoop<Dim> parallelLoop(finePb->_mesh->InteriorFaces);
			parallelLoop.ReserveChunkCoeffsSize(nFaceUnknowns * 2 * nFaceUnknowns);

			parallelLoop.Execute([this, P_algo1, J_faces, nFaceUnknowns](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
				{
					Poisson_HHO_Face<Dim>* fineFace = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
					if (f->IsRemovedOnCoarserGrid)
						chunk->Results.Coeffs.CopyRows(fineFace->Number*nFaceUnknowns, nFaceUnknowns, P_algo1);
					else
						chunk->Results.Coeffs.CopyRows(fineFace->Number*nFaceUnknowns, nFaceUnknowns, J_faces);
				});

			P = SparseMatrix(finePb->HHO->nInteriorFaces * nFaceUnknowns, coarsePb->HHO->nInteriorFaces * nFaceUnknowns);
			parallelLoop.Fill(P);
		}
		else if (_prolongationCode == 3)
		{
			SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
			SparseMatrix J_f_c = GetGlobalCanonicalInjectionMatrixCoarseToFineElements();
			SparseMatrix K_f = GetGlobalMatrixFindFacesWhichReconstructCells(finePb);

			if (ExportMatrices)
			{
				finePb->ExportMatrix(I_c, "I_c");
				finePb->ExportMatrix(J_f_c, "J_f_c");
				finePb->ExportMatrix(K_f, "K_f");
			}

			P = K_f * J_f_c * I_c;
		}
		else if (_prolongationCode == 4)
		{
			int nFaceUnknowns = finePb->HHO->nFaceUnknowns;

			ElementParallelLoop<Dim> parallelLoop(coarsePb->_mesh->Elements);

			parallelLoop.Execute([this, finePb, nFaceUnknowns](Element<Dim>* ce, ParallelChunk<CoeffsChunk>* chunk)
				{
					Poisson_HHO_Element<Dim>* coarseElem = dynamic_cast<Poisson_HHO_Element<Dim>*>(ce);
					DenseMatrix resolveCondensedFinerFacesFromCoarseBoundary = coarseElem->StaticallyCondenseInteriorFinerFaces(this->OperatorMatrix);

					for (int i = 0; i < coarseElem->FinerFacesRemoved.size(); i++)
					{
						Face<Dim>* condensedFineFace = coarseElem->FinerFacesRemoved[i];
						for (Face<Dim>* coarseFace : coarseElem->Faces)
						{
							if (coarseFace->IsDomainBoundary)
								continue;

							BigNumber coarseFaceLocalNumberInCoarseElem = coarseElem->LocalNumberOf(coarseFace);
							chunk->Results.Coeffs.Add(condensedFineFace->Number*nFaceUnknowns, coarseFace->Number*nFaceUnknowns, resolveCondensedFinerFacesFromCoarseBoundary.block(i*nFaceUnknowns, coarseFaceLocalNumberInCoarseElem*nFaceUnknowns, nFaceUnknowns, nFaceUnknowns));
						}
					}

					for (Face<Dim>* cf : coarseElem->Faces)
					{
						if (cf->IsDomainBoundary)
							continue;

						Poisson_HHO_Face<Dim>* coarseFace = dynamic_cast<Poisson_HHO_Face<Dim>*>(cf);
						DenseMatrix local_J_f_c = coarseFace->ComputeCanonicalInjectionMatrixCoarseToFine(finePb->HHO->FaceBasis);
						for (auto fineFace : coarseFace->FinerFaces)
						{
							BigNumber fineFaceLocalNumberInCoarseFace = coarseFace->LocalNumberOf(fineFace);
							chunk->Results.Coeffs.Add(fineFace->Number*nFaceUnknowns, coarseFace->Number*nFaceUnknowns, local_J_f_c.block(fineFaceLocalNumberInCoarseFace*nFaceUnknowns, 0, nFaceUnknowns, nFaceUnknowns));
						}
					}
				});

			P = SparseMatrix(finePb->HHO->nInteriorFaces * nFaceUnknowns, coarsePb->HHO->nInteriorFaces * nFaceUnknowns);
			parallelLoop.Fill(P);
		}
		else if (_prolongationCode == 5)
		{
			SparseMatrix J_faces = GetGlobalCanonicalInjectionMatrixCoarseToFineFaces();
			P = J_faces;
		}
		else
			Utils::FatalError("Unknown prolongationCode.");

		if (ExportMatrices)
			finePb->ExportMatrix(P, "P");
	}

	void SetupRestriction() override
	{
		R = RestrictionScalingFactor() * P.transpose();

		if (ExportMatrices)
			this->_problem->ExportMatrix(R, "R");
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
		/*Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		scalingFactor = coarsePb->_mesh->SkeletonMeasure() / finePb->_mesh->SkeletonMeasure();*/

		//------------- Try 3 -------------//
		// Scaling conservation:
		// Actually rescales correctly, but that's not what we want //
		/*Poisson_HHO<Dim>* finePb = this->_problem;
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
		/*Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		DomFunction constantOne = [](DomPoint p) { return 1; };
		Vector constantOneDoFs = coarsePb->ProjectOnFaceDiscreteSpace(constantOne);

		double coarseProduct = constantOneDoFs.transpose() * coarsePb->A * constantOneDoFs;
		double fineProduct = constantOneDoFs.transpose() * P.transpose() * finePb->A * P * constantOneDoFs;
		scalingFactor = (coarseProduct / fineProduct);

		cout << "scalingFactor = " << scalingFactor << endl;*/
		return scalingFactor;
	}

	inline double Weight(Element<Dim>* element, Face<Dim>* face)
	{
		auto n = element->OuterNormalVector(face);
		double k = (element->DiffTensor * n).dot(n);

		auto n1 = face->Element1->OuterNormalVector(face);
		double k1 = (face->Element1->DiffTensor * n1).dot(n1);

		auto n2 = face->Element2->OuterNormalVector(face);
		double k2 = (face->Element2->DiffTensor * n2).dot(n2);

		return k / (k1 + k2);
	}

	SparseMatrix GetGlobalInterpolationMatrixFromFacesToCells(Poisson_HHO<Dim>* problem)
	{
		int nCellUnknowns = _cellInterpolationBasis->Size();
		int nFaceUnknowns = problem->HHO->nFaceUnknowns;

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nFaceUnknowns, nCellUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);

				DenseMatrix cellInterpMatrix;
				if (_cellInterpolationBasis == _problem->HHO->ReconstructionBasis)
					cellInterpMatrix = element->ReconstructionFromFacesMatrix();
				else if (_cellInterpolationBasis == _problem->HHO->CellBasis)
					cellInterpMatrix = element->SolveCellUnknownsMatrix();
				else
					assert(false);

				for (auto f : element->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

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

	SparseMatrix GetSolveCellUnknownsMatrix(Poisson_HHO<Dim>* problem)
	{
		int nCellUnknowns = problem->HHO->nCellUnknowns;
		int nFaceUnknowns = problem->HHO->nFaceUnknowns;

		//FunctionalBasis<Dim - 1>* faceBasis = problem->HHO->FaceBasis;
		//FunctionalBasis<Dim>* cellInterpolationBasis = _problem->HHO->CellBasis;

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nFaceUnknowns, nCellUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);
				DenseMatrix reconstructMatrix = element->SolveCellUnknownsMatrix();

				for (auto f : element->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

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

	SparseMatrix GetGlobalProjectorMatrixFromCellsToFaces(Poisson_HHO<Dim>* problem)
	{
		int nFaceUnknowns = problem->HHO->nFaceUnknowns;
		int nCellUnknowns = _cellInterpolationBasis->Size();

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nCellUnknowns, nFaceUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);

				for (auto f : element->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

					BigNumber elemGlobalNumber = element->Number;
					BigNumber faceGlobalNumber = face->Number;

					double weight = Weight(element, face);
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, elemGlobalNumber*nCellUnknowns, weight*face->GetProjFromCell(element, _cellInterpolationBasis));
				}
			});

		SparseMatrix Pi(problem->HHO->nTotalFaceUnknowns, problem->HHO->nElements * nCellUnknowns);
		parallelLoop.Fill(Pi);

		return Pi;
	}

	SparseMatrix GetGlobalProjectorMatrixFromCoarseCellsToFineFaces()
	{
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;

		int nFaceUnknowns = finePb->HHO->nFaceUnknowns;
		int nCellUnknowns = _cellInterpolationBasis->Size();

		FaceParallelLoop<Dim> parallelLoop(finePb->_mesh->InteriorFaces);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nCellUnknowns, nFaceUnknowns](Face<Dim>* f, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
				BigNumber faceGlobalNumber = face->Number;
				if (face->IsRemovedOnCoarserGrid)
				{
					Element<Dim>* coarseElem = face->Element1->CoarserElement;
					BigNumber coarseElemGlobalNumber = coarseElem->Number;

					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, coarseElemGlobalNumber*nCellUnknowns, face->GetProjFromCell(coarseElem, _cellInterpolationBasis));
				}
				else
				{
					Element<Dim>* coarseElem1 = face->Element1->CoarserElement;
					Element<Dim>* coarseElem2 = face->Element2->CoarserElement;
					BigNumber coarseElem1GlobalNumber = coarseElem1->Number;
					BigNumber coarseElem2GlobalNumber = coarseElem2->Number;

					double weight1 = Weight(coarseElem1, face->CoarseFace); // Careful with using face->CoarseFace when the mesh isn't nested...
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, coarseElem1GlobalNumber*nCellUnknowns, weight1*face->GetProjFromCell(coarseElem1, _cellInterpolationBasis));

					double weight2 = Weight(coarseElem2, face->CoarseFace);
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, coarseElem2GlobalNumber*nCellUnknowns, weight2*face->GetProjFromCell(coarseElem2, _cellInterpolationBasis));
				}
			});

		SparseMatrix Pi(finePb->HHO->nTotalFaceUnknowns, coarsePb->HHO->nElements * nCellUnknowns);
		parallelLoop.Fill(Pi);

		return Pi;
	}

	// Min 1/2 ||  ||^2
	SparseMatrix GetGlobalMatrixFindFacesWhichReconstructCells(Poisson_HHO<Dim>* problem)
	{
		if (_cellInterpolationBasis != _problem->HHO->ReconstructionBasis)
			assert(false && "This algo is only available if we reconstruct in degree k+1 on the cells.");

		int nFaceUnknowns = problem->HHO->nFaceUnknowns;
		int nCellUnknowns = _cellInterpolationBasis->Size();

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nCellUnknowns, nFaceUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);

				DenseMatrix findFacesMatrix = element->FindFacesPolyWhichReconstructOnTheCell();

				for (auto f : element->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

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
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		Mesh<Dim>* coarseMesh = coarsePb->_mesh;

		int nCellUnknowns = _cellInterpolationBasis->Size();

		ElementParallelLoop<Dim> parallelLoop(coarseMesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nCellUnknowns);

		parallelLoop.Execute([this, nCellUnknowns](Element<Dim>* ce, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* coarseElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(ce);

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

	SparseMatrix GetGlobalCanonicalInjectionMatrixCoarseToFineFaces()
	{
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		Mesh<Dim>* coarseMesh = coarsePb->_mesh;

		FunctionalBasis<Dim - 1>* faceBasis = finePb->HHO->FaceBasis;
		int nFaceUnknowns = faceBasis->Size();

		FaceParallelLoop<Dim> parallelLoop(coarseMesh->InteriorFaces);
		parallelLoop.ReserveChunkCoeffsSize(nFaceUnknowns * 2 * nFaceUnknowns);

		parallelLoop.Execute([this, faceBasis, nFaceUnknowns](Face<Dim>* cf, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Face<Dim>* coarseFace = dynamic_cast<Poisson_HHO_Face<Dim>*>(cf);

				DenseMatrix local_J_f_c = coarseFace->ComputeCanonicalInjectionMatrixCoarseToFine(faceBasis);
				for (auto fineFace : coarseFace->FinerFaces)
				{
					BigNumber coarseFaceGlobalNumber = coarseFace->Number;
					BigNumber fineFaceGlobalNumber = fineFace->Number;
					BigNumber fineFaceLocalNumber = coarseFace->LocalNumberOf(fineFace);

					chunk->Results.Coeffs.Add(fineFaceGlobalNumber*nFaceUnknowns, coarseFaceGlobalNumber*nFaceUnknowns, local_J_f_c.block(fineFaceLocalNumber*nFaceUnknowns, 0, nFaceUnknowns, nFaceUnknowns));
				}
			});

		SparseMatrix J_f_c(finePb->HHO->nInteriorFaces * nFaceUnknowns, coarsePb->HHO->nInteriorFaces * nFaceUnknowns);
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
	Poisson_HHO<Dim>* _problem;
	int _prolongationCode = 1;
	FunctionalBasis<Dim>* _cellInterpolationBasis;
public:

	MultigridForHHO(Poisson_HHO<Dim>* problem, int prolongationCode, FunctionalBasis<Dim>* cellInterpolationBasis)
		: MultigridForHHO(problem, prolongationCode, cellInterpolationBasis, 0)
	{}

	MultigridForHHO(Poisson_HHO<Dim>* problem, int prolongationCode, FunctionalBasis<Dim>* cellInterpolationBasis, int nLevels)
		: Multigrid(nLevels)
	{
		this->_problem = problem;
		this->_prolongationCode = prolongationCode;
		this->_cellInterpolationBasis = cellInterpolationBasis;
		this->BlockSizeForBlockSmoothers = problem->HHO->nFaceUnknowns;
		this->_fineLevel = new LevelForHHO<Dim>(0, problem, _prolongationCode, _cellInterpolationBasis);
	}

	void BeginSerialize(ostream& os) const override
	{
		os << "MultigridForHHO" << endl;
		os << "\t" << "Prolongation       : " << _prolongationCode << endl;
		if (_prolongationCode != 5)
		{
			os << "\t" << "Cell interpolation : ";
			if (_cellInterpolationBasis == _problem->HHO->CellBasis)
				os << "k";
			else if (_cellInterpolationBasis == _problem->HHO->ReconstructionBasis)
				os << "k+1";
			os << endl;
		}
	}
	
protected:
	Level* CreateCoarseLevel(Level* fineLevel) override
	{
		LevelForHHO<Dim>* hhoFineLevel = dynamic_cast<LevelForHHO<Dim>*>(fineLevel);
		Poisson_HHO<Dim>* coarseProblem = hhoFineLevel->_problem->GetProblemOnCoarserMesh();
		LevelForHHO<Dim>* coarseLevel = new LevelForHHO<Dim>(fineLevel->Number + 1, coarseProblem, _prolongationCode, _cellInterpolationBasis);
		return coarseLevel;
	}
};