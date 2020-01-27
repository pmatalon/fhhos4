#pragma once
#include "Multigrid.h"
#include "SmootherFactory.h"
#include "../HHO/Poisson_HHO.h"
#include "../Utils/ElementParallelLoop.h"
using namespace std;

template <int Dim>
class LevelForHHO : public Level
{
private:
	Poisson_HHO<Dim>* _problem;
	int _algoNumber = 1;
	FunctionalBasis<Dim>* _cellInterpolationBasis;
public:
	bool ExportMatrices = false;

	LevelForHHO(int number, Poisson_HHO<Dim>* problem, bool useGalerkinOperator, int algoNumber, FunctionalBasis<Dim>* cellInterpolationBasis)
		: Level(number)
	{
		this->UseGalerkinOperator = useGalerkinOperator;
		this->_problem = problem;
		this->_algoNumber = algoNumber;
		this->_cellInterpolationBasis = cellInterpolationBasis;

		//problem->ExportFaces("L" + to_string(number));
	}

	void ExportVector(Vector& v, string suffix)
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
		cout << "\t\tMesh              : " << this->_problem->_mesh->Elements.size() << " elements" << endl;
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

		SparseMatrix I_c = GetGlobalInterpolationMatrixFromFacesToCells(coarsePb);
		SparseMatrix J_f_c = GetGlobalCanonicalInjectionMatrixCoarseToFineElements();
		SparseMatrix Pi_f = GetGlobalProjectorMatrixFromCellsToFaces(finePb);

		if (ExportMatrices)
		{
			SparseMatrix SolveCellUnknowns = GetSolveCellUnknownsMatrix(coarsePb);
			finePb->ExportMatrix(I_c, "I_c");
			finePb->ExportMatrix(J_f_c, "J_f_c");
			finePb->ExportMatrix(Pi_f, "Pi_f");
			finePb->ExportMatrix(SolveCellUnknowns, "SolveCellUnknowns");
		}

		SparseMatrix P_algo1 = Pi_f * J_f_c * I_c;

		if (_algoNumber == 1)
			P = P_algo1;
		else
		{
			int nFaceUnknowns = finePb->HHO->nFaceUnknowns;
			SparseMatrix J_faces = GetGlobalCanonicalInjectionMatrixCoarseToFineFaces();

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
				for (auto f : element->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

					BigNumber elemGlobalNumber = element->Number;
					BigNumber faceGlobalNumber = face->Number;
					BigNumber faceLocalNumber = element->LocalNumberOf(face);

					DenseMatrix cellInterpMatrix;
					if (_cellInterpolationBasis == _problem->HHO->ReconstructionBasis)
						cellInterpMatrix = element->ReconstructionFromFacesMatrix();
					else if (_cellInterpolationBasis == _problem->HHO->CellBasis)
						cellInterpMatrix = element->SolveCellUnknownsMatrix();
					else
						assert(false);

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

		FunctionalBasis<Dim - 1>* faceBasis = problem->HHO->FaceBasis;
		FunctionalBasis<Dim>* cellInterpolationBasis = _problem->HHO->CellBasis;

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nFaceUnknowns, nCellUnknowns](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* element = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);
				for (auto f : element->Faces)
				{
					if (f->IsDomainBoundary)
						continue;

					Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

					BigNumber elemGlobalNumber = element->Number;
					BigNumber faceGlobalNumber = face->Number;
					BigNumber faceLocalNumber = element->LocalNumberOf(face);
					DenseMatrix reconstructMatrix = element->SolveCellUnknownsMatrix();
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
					BigNumber faceLocalNumber = element->LocalNumberOf(face);

					double weight = Weight(element, face);
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, elemGlobalNumber*nCellUnknowns, weight*face->GetProjFromCell(element, _cellInterpolationBasis));
				}
			});

		SparseMatrix Pi(problem->HHO->nTotalFaceUnknowns, problem->HHO->nElements * nCellUnknowns);
		parallelLoop.Fill(Pi);

		return Pi;
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
	bool _automaticNumberOfLevels;
	int _nLevels;
	int _algoNumber = 1;
	FunctionalBasis<Dim>* _cellInterpolationBasis;
public:
	int MatrixMaxSizeForCoarsestLevel;
	bool UseGalerkinOperator = false;
	string PreSmootherCode = "bgs";
	string PostSmootherCode = "rbgs";
	int PreSmoothingIterations = 1;
	int PostSmoothingIterations = 1;
	CoarseningStrategy CoarseningStgy = CoarseningStrategy::Standard;
	bool ExportMatrices = false;

	MultigridForHHO(Poisson_HHO<Dim>* problem, int algoNumber, FunctionalBasis<Dim>* cellInterpolationBasis) : MultigridForHHO(problem, algoNumber, cellInterpolationBasis, 0)
	{}

	MultigridForHHO(Poisson_HHO<Dim>* problem, int algoNumber, FunctionalBasis<Dim>* cellInterpolationBasis, int nLevels) : Multigrid()
	{
		this->_problem = problem;
		this->_algoNumber = algoNumber;
		this->_cellInterpolationBasis = cellInterpolationBasis;
		this->_automaticNumberOfLevels = nLevels <= 0;
		
		if (this->_automaticNumberOfLevels)
		{
			this->_nLevels = 0;
			this->MatrixMaxSizeForCoarsestLevel = 1000;
		}
		else
		{
			this->_nLevels = nLevels;
			this->MatrixMaxSizeForCoarsestLevel = 0;
		}

		this->_fineLevel = new LevelForHHO<Dim>(0, problem, false, _algoNumber, _cellInterpolationBasis);
	}
	
	void Serialize(ostream& os) const override
	{
		os << "MultigridForHHO" << endl;
		os << "\t" << "Cycle              : ";
		if (this->WLoops == 1)
			os << "V-cycle" << endl;
		else
			os << "W-cycle (" << this->WLoops << " loops)" << endl;
		os << "\t" << "Levels             : ";
		if (_automaticNumberOfLevels && _nLevels == 0)
			os << "automatic coarsening until matrix size <= " << MatrixMaxSizeForCoarsestLevel << endl;
		else if (_automaticNumberOfLevels && _nLevels > 0)
			os << _nLevels << " (automatic)" << endl;
		else
			os << _nLevels << endl;
		os << "\t" << "Coarsening strategy: ";
		if (CoarseningStgy == CoarseningStrategy::Standard)
			os << "standard" << endl;
		else if (CoarseningStgy == CoarseningStrategy::Agglomeration)
			os << "agglomeration" << endl;
		else if (CoarseningStgy == CoarseningStrategy::StructuredRefinement)
			os << "structured refinement from coarse mesh" << endl;
		Smoother* preSmoother = SmootherFactory::Create(PreSmootherCode, PreSmoothingIterations, _problem->HHO->nFaceUnknowns);
		Smoother* postSmoother = SmootherFactory::Create(PostSmootherCode, PostSmoothingIterations, _problem->HHO->nFaceUnknowns);
		os << "\t" << "Pre-smoothing      : " << *preSmoother << endl;
		os << "\t" << "Post-smoothing     : " << *postSmoother;
		delete preSmoother;
		delete postSmoother;
	}

	void Setup(const SparseMatrix& A) override
	{
		IterativeSolver::Setup(A);

		cout << "Setup..." << endl;

		LevelForHHO<Dim>* finerLevel = dynamic_cast<LevelForHHO<Dim>*>(this->_fineLevel);
		finerLevel->OperatorMatrix = A;
		finerLevel->PreSmoother = SmootherFactory::Create(PreSmootherCode, PreSmoothingIterations, _problem->HHO->nFaceUnknowns);
		finerLevel->PostSmoother = SmootherFactory::Create(PostSmootherCode, PostSmoothingIterations, _problem->HHO->nFaceUnknowns);
		finerLevel->ExportMatrices = ExportMatrices;
		Poisson_HHO<Dim>* problem = _problem;

		bool noCoarserMeshProvided = false;
		bool coarsestPossibleMeshReached = false;

		int levelNumber = 0;
		while ((_automaticNumberOfLevels && problem->HHO->nTotalFaceUnknowns > MatrixMaxSizeForCoarsestLevel) || (levelNumber < _nLevels - 1))
		{
			// Can we coarsen the mesh?
			if (this->CoarseningStgy == CoarseningStrategy::StructuredRefinement && problem->_mesh->CoarseMesh == nullptr)
			{
				noCoarserMeshProvided = true;
				break;
			}
			problem->_mesh->CoarsenMesh(this->CoarseningStgy);
			if (problem->_mesh->CoarseMesh->InteriorFaces.size() == 0)
			{
				coarsestPossibleMeshReached = true;
				break;
			}

			// Build coarse level
			problem = problem->GetProblemOnCoarserMesh();
			levelNumber++;
			LevelForHHO<Dim>* coarseLevel = new LevelForHHO<Dim>(levelNumber, problem, UseGalerkinOperator, _algoNumber, _cellInterpolationBasis);
			coarseLevel->PreSmoother = SmootherFactory::Create(PreSmootherCode, PreSmoothingIterations, problem->HHO->nFaceUnknowns);
			coarseLevel->PostSmoother = SmootherFactory::Create(PostSmootherCode, PostSmoothingIterations, problem->HHO->nFaceUnknowns);
			coarseLevel->ExportMatrices = ExportMatrices;

			// Link between levels
			finerLevel->CoarserLevel = coarseLevel;
			coarseLevel->FinerLevel = finerLevel;

			// Setup fine level
			finerLevel->Setup();

			finerLevel = coarseLevel;
		}

		_nLevels = levelNumber + 1;
		finerLevel->Setup();

		if (coarsestPossibleMeshReached)
			cout << Utils::BeginYellow << "Warning: impossible to coarsen the mesh any more." << Utils::EndColor;
		if (noCoarserMeshProvided)
			cout << Utils::BeginYellow << "Warning: cannot build coarser level because no coarser mesh has been provided." << Utils::EndColor;

		this->SetupCoarseSolver();
		cout << "\t--> " << _nLevels << " levels built." << endl;
		
		if (this->WLoops > 1)
			PrintCycleSchema();
	}
};