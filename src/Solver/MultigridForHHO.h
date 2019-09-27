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
public:
	LevelForHHO(int number, Poisson_HHO<Dim>* problem, bool useGalerkinOperator, int algoNumber)
		: Level(number)
	{
		this->UseGalerkinOperator = useGalerkinOperator;
		this->_problem = problem;
		this->_algoNumber = algoNumber;

		//problem->ExportFaces("L" + to_string(number));
	}

	void ExportVector(Eigen::VectorXd& v, string suffix)
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

		/*finePb->ExportMatrix(I_c, "I_c");
		finePb->ExportMatrix(J_f_c, "J_f_c");
		finePb->ExportMatrix(Pi_f, "Pi_f");*/

		SparseMatrix P_algo1 = Pi_f * J_f_c * I_c;

		if (_algoNumber == 1)
			P = P_algo1;
		else
		{
			int nFaceUnknowns = finePb->HHO.nLocalFaceUnknowns;
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

			P = SparseMatrix(finePb->HHO.nInteriorFaces * nFaceUnknowns, coarsePb->HHO.nInteriorFaces * nFaceUnknowns);
			parallelLoop.Fill(P);
		}
		//finePb->ExportMatrix(P, "P");
	}

	void SetupRestriction() override
	{
		R = P.transpose();
		//this->_problem->ExportMatrix(R, "R");
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

	FunctionalBasis<Dim>* GetCellInterpolationBasis(Poisson_HHO<Dim>* problem)
	{
		FunctionalBasis<Dim-1>* faceBasis = dynamic_cast<Poisson_HHO_Element<Dim>*>(problem->_mesh->Elements[0])->FaceBasis;
		int faceDegreePSpace = faceBasis->GetDegree();

		//int cellDegreeQSpace = (int)floor(2 * sqrt(faceDegreePSpace + 1) - 1);
		//FunctionalBasis<Dim>* cellInterpolationBasis = new FunctionalBasis<Dim>("monomials", cellDegreeQSpace, true);
		int cellDegreePSpace = faceDegreePSpace + 1;
		FunctionalBasis<Dim>* cellInterpolationBasis = new FunctionalBasis<Dim>(faceBasis->BasisCode(), cellDegreePSpace, false);

		return cellInterpolationBasis;
	}

	SparseMatrix GetGlobalInterpolationMatrixFromFacesToCells(Poisson_HHO<Dim>* problem)
	{
		return GetReconstructionMatrix(problem);
	}

	SparseMatrix GetReconstructionMatrix(Poisson_HHO<Dim>* problem)
	{
		int nCellUnknowns = problem->HHO.nLocalCellUnknowns;
		int nFaceUnknowns = problem->HHO.nLocalFaceUnknowns;

		FunctionalBasis<Dim>* cellInterpolationBasis = GetCellInterpolationBasis(problem);

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nFaceUnknowns, cellInterpolationBasis](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
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
					Eigen::MatrixXd reconstructMatrix = element->ReconstructionFromFacesMatrix();
					chunk->Results.Coeffs.Add(elemGlobalNumber * cellInterpolationBasis->Size(), faceGlobalNumber * nFaceUnknowns, reconstructMatrix.block(0, faceLocalNumber*nFaceUnknowns, cellInterpolationBasis->Size(), nFaceUnknowns));
				}
			});

		SparseMatrix M(problem->HHO.nElements * cellInterpolationBasis->Size(), problem->HHO.nTotalFaceUnknowns);
		parallelLoop.Fill(M);
		delete cellInterpolationBasis;

		return M;
	}

	SparseMatrix GetGlobalProjectorMatrixFromCellsToFaces(Poisson_HHO<Dim>* problem)
	{
		int nFaceUnknowns = problem->HHO.nLocalFaceUnknowns;

		FunctionalBasis<Dim>* cellInterpolationBasis = GetCellInterpolationBasis(problem);
		int nCellUnknowns = cellInterpolationBasis->Size();

		ElementParallelLoop<Dim> parallelLoop(problem->_mesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nFaceUnknowns);

		parallelLoop.Execute([this, nCellUnknowns, nFaceUnknowns, cellInterpolationBasis](Element<Dim>* e, ParallelChunk<CoeffsChunk>* chunk)
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
					chunk->Results.Coeffs.Add(faceGlobalNumber*nFaceUnknowns, elemGlobalNumber*nCellUnknowns, weight*face->GetProjFromCell(element, cellInterpolationBasis));
				}
			});

		SparseMatrix Pi(problem->HHO.nTotalFaceUnknowns, problem->HHO.nElements * nCellUnknowns);
		parallelLoop.Fill(Pi);

		return Pi;
	}

	SparseMatrix GetGlobalCanonicalInjectionMatrixCoarseToFineElements()
	{
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		Mesh<Dim>* coarseMesh = coarsePb->_mesh;

		FunctionalBasis<Dim>* cellInterpolationBasis = GetCellInterpolationBasis(finePb);
		int nCellUnknowns = cellInterpolationBasis->Size();

		ElementParallelLoop<Dim> parallelLoop(coarseMesh->Elements);
		parallelLoop.ReserveChunkCoeffsSize(nCellUnknowns * 4 * nCellUnknowns);

		parallelLoop.Execute([this, cellInterpolationBasis, nCellUnknowns](Element<Dim>* ce, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Element<Dim>* coarseElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(ce);

				Eigen::MatrixXd local_J_f_c = coarseElement->ComputeCanonicalInjectionMatrixCoarseToFine(cellInterpolationBasis);
				for (auto fineElement : coarseElement->FinerElements)
				{
					BigNumber coarseElemGlobalNumber = coarseElement->Number;
					BigNumber fineElemGlobalNumber = fineElement->Number;
					BigNumber fineElemLocalNumber = coarseElement->LocalNumberOf(fineElement);

					chunk->Results.Coeffs.Add(fineElemGlobalNumber*nCellUnknowns, coarseElemGlobalNumber*nCellUnknowns, local_J_f_c.block(fineElemLocalNumber*nCellUnknowns, 0, nCellUnknowns, nCellUnknowns));
				}
			});

		SparseMatrix J_f_c(finePb->HHO.nElements * nCellUnknowns, coarsePb->HHO.nElements * nCellUnknowns);
		parallelLoop.Fill(J_f_c);
		return J_f_c;
	}

	SparseMatrix GetGlobalCanonicalInjectionMatrixCoarseToFineFaces()
	{
		Poisson_HHO<Dim>* finePb = this->_problem;
		Poisson_HHO<Dim>* coarsePb = dynamic_cast<LevelForHHO<Dim>*>(CoarserLevel)->_problem;
		Mesh<Dim>* coarseMesh = coarsePb->_mesh;

		FunctionalBasis<Dim - 1>* faceBasis = dynamic_cast<Poisson_HHO_Face<Dim>*>(finePb->_mesh->Faces[0])->FaceBasis;
		int nFaceUnknowns = faceBasis->Size();

		FaceParallelLoop<Dim> parallelLoop(coarseMesh->InteriorFaces);
		parallelLoop.ReserveChunkCoeffsSize(nFaceUnknowns * 2 * nFaceUnknowns);

		parallelLoop.Execute([this, faceBasis, nFaceUnknowns](Face<Dim>* cf, ParallelChunk<CoeffsChunk>* chunk)
			{
				Poisson_HHO_Face<Dim>* coarseFace = dynamic_cast<Poisson_HHO_Face<Dim>*>(cf);

				Eigen::MatrixXd local_J_f_c = coarseFace->ComputeCanonicalInjectionMatrixCoarseToFine(faceBasis);
				for (auto fineFace : coarseFace->FinerFaces)
				{
					BigNumber coarseFaceGlobalNumber = coarseFace->Number;
					BigNumber fineFaceGlobalNumber = fineFace->Number;
					BigNumber fineFaceLocalNumber = coarseFace->LocalNumberOf(fineFace);

					chunk->Results.Coeffs.Add(fineFaceGlobalNumber*nFaceUnknowns, coarseFaceGlobalNumber*nFaceUnknowns, local_J_f_c.block(fineFaceLocalNumber*nFaceUnknowns, 0, nFaceUnknowns, nFaceUnknowns));
				}
			});

		SparseMatrix J_f_c(finePb->HHO.nInteriorFaces * nFaceUnknowns, coarsePb->HHO.nInteriorFaces * nFaceUnknowns);
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
public:
	int MatrixMaxSizeForCoarsestLevel;
	bool UseGalerkinOperator = false;
	string PreSmootherCode = "bgs";
	string PostSmootherCode = "rbgs";
	int PreSmoothingIterations = 1;
	int PostSmoothingIterations = 1;
	CoarseningStrategy CoarseningStgy = CoarseningStrategy::Standard;

	MultigridForHHO(Poisson_HHO<Dim>* problem, int algoNumber) : MultigridForHHO(problem, algoNumber, 0)
	{}

	MultigridForHHO(Poisson_HHO<Dim>* problem, int algoNumber, int nLevels) : Multigrid()
	{
		this->_problem = problem;
		this->_algoNumber = algoNumber;
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

		this->_fineLevel = new LevelForHHO<Dim>(0, problem, false, _algoNumber);
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
		Smoother* preSmoother = SmootherFactory::Create(PreSmootherCode, PreSmoothingIterations, _problem->HHO.nLocalFaceUnknowns);
		Smoother* postSmoother = SmootherFactory::Create(PostSmootherCode, PostSmoothingIterations, _problem->HHO.nLocalFaceUnknowns);
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
		finerLevel->PreSmoother = SmootherFactory::Create(PreSmootherCode, PreSmoothingIterations, _problem->HHO.nLocalFaceUnknowns);
		finerLevel->PostSmoother = SmootherFactory::Create(PostSmootherCode, PostSmoothingIterations, _problem->HHO.nLocalFaceUnknowns);
		Poisson_HHO<Dim>* problem = _problem;
		if (!_automaticNumberOfLevels)
		{
			for (int levelNumber = 1; levelNumber < _nLevels; levelNumber++)
			{
				//cout << "fine mesh" << endl << *(problem->_mesh) << endl << endl;
				problem->_mesh->CoarsenMesh(this->CoarseningStgy);
				if (problem->_mesh->CoarseMesh->InteriorFaces.size() == 0)
				{
					cout << "Warning: impossible to coarsen the mesh any more.";
					break;
				}
				problem = problem->GetProblemOnCoarserMesh();
				//cout << "coarse mesh" << endl << *(problem->_mesh) << endl << endl;
				LevelForHHO<Dim>* coarseLevel = new LevelForHHO<Dim>(levelNumber, problem, UseGalerkinOperator, _algoNumber);
				coarseLevel->PreSmoother = SmootherFactory::Create(PreSmootherCode, PreSmoothingIterations, problem->HHO.nLocalFaceUnknowns);
				coarseLevel->PostSmoother = SmootherFactory::Create(PostSmootherCode, PostSmoothingIterations, problem->HHO.nLocalFaceUnknowns);

				finerLevel->CoarserLevel = coarseLevel;
				coarseLevel->FinerLevel = finerLevel;

				finerLevel->Setup();

				finerLevel = coarseLevel;
			}
			finerLevel->Setup();
			this->SetupCoarseSolver();
		}
		else
		{
			int levelNumber = 0;
			while (problem->HHO.nTotalFaceUnknowns > MatrixMaxSizeForCoarsestLevel)
			{
				problem->_mesh->CoarsenMesh(this->CoarseningStgy);
				auto nCoarseInteriorFaces = problem->_mesh->CoarseMesh->InteriorFaces.size();
				auto nCoarseUnknowns = nCoarseInteriorFaces * problem->HHO.nLocalFaceUnknowns;
				if (nCoarseInteriorFaces == 0)
					break;
				problem = problem->GetProblemOnCoarserMesh();
				levelNumber++;
				LevelForHHO<Dim>* coarseLevel = new LevelForHHO<Dim>(levelNumber, problem, UseGalerkinOperator, _algoNumber);
				coarseLevel->PreSmoother = SmootherFactory::Create(PreSmootherCode, PreSmoothingIterations, problem->HHO.nLocalFaceUnknowns);
				coarseLevel->PostSmoother = SmootherFactory::Create(PostSmootherCode, PostSmoothingIterations, problem->HHO.nLocalFaceUnknowns);

				finerLevel->CoarserLevel = coarseLevel;
				coarseLevel->FinerLevel = finerLevel;

				finerLevel->Setup();

				finerLevel = coarseLevel;
			}
			finerLevel->Setup();
			_nLevels = levelNumber + 1;
			this->SetupCoarseSolver();
			cout << "\t--> " << _nLevels << " levels built." << endl;
		}
		
		if (this->WLoops > 1)
			PrintCycleSchema();
	}
};