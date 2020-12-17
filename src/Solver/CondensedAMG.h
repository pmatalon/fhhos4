#pragma once
#include <mutex>
#include "CondensedAlgebraicMesh.h"
#include "Multigrid.h"
using namespace std;

class CondensedLevel : public Level
{
private:
	CAMGProlongation _prolongation = CAMGProlongation::P;
	int _cellBlockSize;
	int _faceBlockSize;
	CondensedAlgebraicMesh _mesh;
public:
	const SparseMatrix* A_T_T;
	const SparseMatrix* A_T_F;
	const SparseMatrix* A_F_F;
	const SparseMatrix* inv_A_T_T;
private:
	SparseMatrix Q_T;
	SparseMatrix Q_F;

	SparseMatrix A_T_Tc;
	SparseMatrix A_T_Fc;
	SparseMatrix A_F_Fc;
	SparseMatrix* inv_A_T_Tc;

public:
	CondensedLevel(int number, int cellBlockSize, int faceBlockSize, CAMGProlongation prolongation)
		: Level(number), _mesh(cellBlockSize, faceBlockSize)
	{
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
		this->_prolongation = prolongation;
	}

	BigNumber NUnknowns() override
	{
		return A_F_F->rows();
	}

	void ExportVector(const Vector& v, string suffix, int levelNumber) override
	{ }
	void ExportMatrix(const SparseMatrix& M, string suffix, int levelNumber) override
	{ }

	void CoarsenMesh(CoarseningStrategy coarseningStgy, int coarseningFactor, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) override
	{
		cout << "\tDouble pairwise aggregation" << endl;

		_mesh.Build(*A_T_T, *A_T_F);

		//----------------------------//
		// First pairwise aggregation //
		//----------------------------//

		_mesh.PairWiseAggregate(coarsestPossibleMeshReached);
		if (coarsestPossibleMeshReached)
			return;

		/*cout << "------------- Elem" << endl;
		for (ElementAggregate& agg : _mesh._coarseElements)
		{
			cout << "(";
			for (auto e : agg.FineElements)
				cout << e->Number << ", ";
			cout << ")" << endl;
		}

		cout << "------------- Faces" << endl;
		for (FaceAggregate& agg : _mesh._coarseFaces)
		{
			cout << "(";
			for (auto f : agg.FineFaces)
				cout << f->Number << ", ";
			cout << ")" << endl;
		}*/

		// Cell-prolongation operator with only one 1 coefficient per row
		SparseMatrix Q_T1 = BuildQ_T(_mesh);
		// Cell-prolongation operator with only one 1 coefficient per row for kept or aggregated faces, average for removed faces
		SparseMatrix Q_F1 = BuildQ_F(_mesh, *A_F_F);

		// Intermediate coarse operators
		SparseMatrix A_T_T1 = Q_T1.transpose() * (*A_T_T) * Q_T1;
		SparseMatrix A_T_F1 = Q_T1.transpose() * (*A_T_F) * Q_F1;
		SparseMatrix A_F_F1 = Q_F1.transpose() * (*A_F_F) * Q_F1;

		//-----------------------------//
		// Second pairwise aggregation //
		//-----------------------------//

		CondensedAlgebraicMesh coarseMesh(_cellBlockSize, _faceBlockSize);
		coarseMesh.Build(A_T_T1, A_T_F1);
		coarseMesh.PairWiseAggregate(coarsestPossibleMeshReached);
		if (coarsestPossibleMeshReached)
			return;

		SparseMatrix Q_T2 = BuildQ_T(coarseMesh);
		SparseMatrix Q_F2 = BuildQ_F(coarseMesh, A_F_F1);

		this->A_T_Tc = Q_T2.transpose() * A_T_T1 * Q_T2;
		this->A_T_Fc = Q_T2.transpose() * A_T_F1 * Q_F2;

		/*
		// Assembly of the double-aggregates
		vector<vector<AlgebraicElement*>> doubleElemAggregates(coarseMesh._coarseElements.size());
		NumberParallelLoop<EmptyResultChunk> parallelLoopE(coarseMesh._coarseElements.size());
		parallelLoopE.Execute([this, &coarseMesh, &doubleElemAggregates](BigNumber secondAggregNumber)
			{
				ElementAggregate& agg2 = coarseMesh._coarseElements[secondAggregNumber];
				vector<AlgebraicElement*> doubleAggregate;
				for (AlgebraicElement* ce1 : agg2.FineElements)
				{
					ElementAggregate& agg1 = _mesh._coarseElements[ce1->Number];
					for (AlgebraicElement* fe : agg1.FineElements)
						doubleAggregate.push_back(fe);
				}
				doubleElemAggregates[secondAggregNumber] = doubleAggregate;
			});

		vector<vector<AlgebraicFace*>> doubleFaceAggregates(coarseMesh._coarseFaces.size());
		NumberParallelLoop<EmptyResultChunk> parallelLoopF(coarseMesh._coarseFaces.size());
		parallelLoopF.Execute([this, &coarseMesh, &doubleFaceAggregates](BigNumber secondAggregNumber)
			{
				FaceAggregate& agg2 = coarseMesh._coarseFaces[secondAggregNumber];
				vector<AlgebraicFace*> doubleAggregate;
				for (AlgebraicFace* cf1 : agg2.FineFaces)
				{
					FaceAggregate& agg1 = _mesh._coarseFaces[cf1->Number];
					for (AlgebraicFace* fe : agg1.FineFaces)
						doubleAggregate.push_back(fe);
				}
				doubleFaceAggregates[secondAggregNumber] = doubleAggregate;
			});
		
		SparseMatrix Q_T = BuildQ_T(doubleElemAggregates);
		SparseMatrix Q_F = BuildQ_F(doubleFaceAggregates);
		*/

		//------------------//
		// Coarse operators //
		//------------------//

		SparseMatrix Q_T = Q_T1 * Q_T2;
		SparseMatrix Q_F = Q_F1 * Q_F2;

		//-----------------------//
		// Prolongation operator //
		//-----------------------//

		if (_prolongation == CAMGProlongation::Q_F)
		{
			this->P = Q_F;
		}
		else if (_prolongation == CAMGProlongation::P)
		{
			SparseMatrix inv_A_T_T1 = Utils::InvertBlockDiagMatrix(A_T_T1, _cellBlockSize);
			// Theta: reconstruction from the coarse faces to the coarse cells
			SparseMatrix Theta1 = -inv_A_T_T1 * A_T_F1;
			// Pi: average on both sides of each face
			SparseMatrix Pi1 = BuildTrace(_mesh);

			this->inv_A_T_Tc = new SparseMatrix(Utils::InvertBlockDiagMatrix(A_T_Tc, _cellBlockSize));
			SparseMatrix Theta2 = -(*inv_A_T_Tc) * A_T_Fc;
			SparseMatrix Pi2 = BuildTrace(coarseMesh);

			this->P = Pi2 * Q_T * Theta2;
		}
		else if (_prolongation == CAMGProlongation::P1P2)
		{
			SparseMatrix inv_A_T_T1 = Utils::InvertBlockDiagMatrix(A_T_T1, _cellBlockSize);
			// Theta: reconstruction from the coarse faces to the coarse cells
			SparseMatrix Theta1 = -inv_A_T_T1 * A_T_F1;
			// Pi: average on both sides of each face
			SparseMatrix Pi1 = BuildTrace(_mesh);
			SparseMatrix P1 = Pi1 * Q_T1 * Theta1;

			this->inv_A_T_Tc = new SparseMatrix(Utils::InvertBlockDiagMatrix(A_T_Tc, _cellBlockSize));
			SparseMatrix Theta2 = -(*inv_A_T_Tc) * A_T_Fc;
			SparseMatrix Pi2 = BuildTrace(coarseMesh);

			SparseMatrix P2 = Pi2 * Q_T2 * Theta2;
			this->P = P1 * P2;
		}
		else
			Utils::FatalError("Unmanaged prolongation");

		cout << "Computing coarse A_F_F" << endl;

		//if (!this->UseGalerkinOperator)
		//{
			this->A_F_Fc = Q_F.transpose() * (*A_F_F) * Q_F;
			//this->A_F_Fc = P.transpose() * (*A_F_F) * P;
		//}
	}

private:
	// Cell prolongation Q_T with only one 1 coefficient per row
	SparseMatrix BuildQ_T(const CondensedAlgebraicMesh& mesh)
	{
		DenseMatrix Id = DenseMatrix::Identity(_cellBlockSize, _cellBlockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoopQ_T(mesh._elements.size());
		parallelLoopQ_T.Execute([this, &mesh, &Id](BigNumber elemNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const AlgebraicElement& elem = mesh._elements[elemNumber];
				chunk->Results.Coeffs.Add(elem.Number, elem.CoarseElement->Number, Id);
			});
		SparseMatrix Q_T = SparseMatrix(mesh._elements.size()*_cellBlockSize, mesh._coarseElements.size()*_cellBlockSize);
		parallelLoopQ_T.Fill(Q_T);
		return Q_T;
	}

	// Face prolongation Q_F
	SparseMatrix BuildQ_F(const CondensedAlgebraicMesh& mesh, const SparseMatrix A_F_F)
	{
		bool enableAnisotropyManagement = false;
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);

		NumberParallelLoop<CoeffsChunk> parallelLoop1(mesh._faces.size());
		parallelLoop1.Execute([this, &mesh, &Id, enableAnisotropyManagement, &A_F_F](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const AlgebraicFace* face = &mesh._faces[faceNumber];
				if (face->IsRemovedOnCoarseMesh)
				{
					ElementAggregate* elemAggreg = face->Elements[0]->CoarseElement;

					if (enableAnisotropyManagement)
					{
						map<FaceAggregate*, double> couplings;
						double totalCouplings = 0;
						for (FaceAggregate* coarseFace : elemAggreg->CoarseFaces)
						{
							double avgCouplingCoarseFace = 0;
							for (AlgebraicFace* f : coarseFace->FineFaces)
							{
								DenseMatrix couplingBlock = A_F_F.block(face->Number*_faceBlockSize, f->Number*_faceBlockSize, _faceBlockSize, _faceBlockSize);
								double couplingFineFace = couplingBlock(0, 0);
								avgCouplingCoarseFace += couplingFineFace;
							}
							avgCouplingCoarseFace /= coarseFace->FineFaces.size();
							if (avgCouplingCoarseFace < 0)
							{
								couplings.insert({ coarseFace, avgCouplingCoarseFace });
								totalCouplings += avgCouplingCoarseFace;
							}
						}

						for (auto it = couplings.begin(); it != couplings.end(); it++)
						{
							FaceAggregate* coarseFace = it->first;
							double coupling = it->second;
							chunk->Results.Coeffs.Add(face->Number, coarseFace->Number, -coupling / abs(totalCouplings) * Id);
						}
					}
					else
					{
						for (FaceAggregate* coarseFace : elemAggreg->CoarseFaces)
							chunk->Results.Coeffs.Add(face->Number, coarseFace->Number, 1.0 / elemAggreg->CoarseFaces.size() * Id);
					}
				}
				else
					chunk->Results.Coeffs.Add(face->Number, face->CoarseFace->Number, Id);
			});


		SparseMatrix Q_F = SparseMatrix(mesh._faces.size()*_faceBlockSize, mesh._coarseFaces.size()*_faceBlockSize);
		parallelLoop1.Fill(Q_F);
		return Q_F;
	}

	SparseMatrix BuildTrace(const CondensedAlgebraicMesh& mesh)
	{
		// Pi: average on both sides of each face
		DenseMatrix traceOfConstant = DenseMatrix::Zero(_faceBlockSize, _cellBlockSize);
		traceOfConstant(0, 0) = 1;

		NumberParallelLoop<CoeffsChunk> parallelLoopPi(mesh._faces.size());
		parallelLoopPi.Execute([&mesh, &traceOfConstant](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const AlgebraicFace& face = mesh._faces[faceNumber];
				if (face.IsRemovedOnCoarseMesh)
					chunk->Results.Coeffs.Add(faceNumber, face.Elements[0]->Number, traceOfConstant);
				else
				{
					for (AlgebraicElement* elem : face.Elements)
						chunk->Results.Coeffs.Add(faceNumber, elem->Number, 1.0 / face.Elements.size()*traceOfConstant);
				}
			});
		SparseMatrix Pi = SparseMatrix(mesh._faces.size()*_faceBlockSize, mesh._elements.size()*_cellBlockSize);
		parallelLoopPi.Fill(Pi);
		return Pi;
	}

/*
	// Cell prolongation Q_T with only one 1 coefficient per row
	SparseMatrix BuildQ_T(const vector<vector<AlgebraicElement*>>& aggregates)
	{
		DenseMatrix Id = DenseMatrix::Identity(_cellBlockSize, _cellBlockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoop(aggregates.size());
		parallelLoop.Execute([this, &aggregates, &Id](BigNumber aggregNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const vector<AlgebraicElement*>& aggreg = aggregates[aggregNumber];
				for (AlgebraicElement* e : aggreg)
					chunk->Results.Coeffs.Add(e->Number, aggregNumber, Id);
			});
		SparseMatrix Q_T = SparseMatrix(_mesh._elements.size()*_cellBlockSize, aggregates.size()*_cellBlockSize);
		parallelLoop.Fill(Q_T);
		return Q_T;
	}

	// Face prolongation Q_F with only one 1 coefficient per row
	SparseMatrix BuildQ_F(vector<vector<AlgebraicFace*>> aggregates)
	{
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoop(_mesh._faces.size());
		parallelLoop.Execute([this, &Id](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const AlgebraicFace* face = &_mesh._faces[faceNumber];
				if (face->IsRemovedOnCoarseMesh)
				{
					ElementAggregate* elemAggreg = face->Elements[0]->CoarseElement;
					double sum = 0;
					for (FaceAggregate* coarseFace : elemAggreg->CoarseFaces)
						chunk->Results.Coeffs.Add(face->Number, coarseFace->Number, 1.0 / elemAggreg->CoarseFaces.size() * Id);
				}
				else
					chunk->Results.Coeffs.Add(face->Number, face->CoarseFace->Number, Id);
			});

		SparseMatrix Q_F = SparseMatrix(_mesh._faces.size()*_faceBlockSize, aggregates.size()*_faceBlockSize);
		parallelLoop.Fill(Q_F);
		return Q_F;
	}*/
public:
	void SetupDiscretizedOperator() override 
	{
		SparseMatrix* schur = new SparseMatrix(*A_F_F - (A_T_F->transpose()) * (*inv_A_T_T) * (*A_T_F));
		this->OperatorMatrix = schur;
	}

	void OnStartSetup() override
	{
		cout << "\t\tMesh                : " << this->A_T_T->rows() / _cellBlockSize << " elements, " << this->A_T_F->cols() / _faceBlockSize << " faces";
		if (!this->IsFinestLevel())
		{
			CondensedLevel* fine = dynamic_cast<CondensedLevel*>(this->FinerLevel);
			double nFine = fine->A_T_F->cols();
			double nCoarse = this->A_T_F->cols();
			cout << ", coarsening factor = " << (nFine/nCoarse);
		}
		cout << endl;
	}

	void SetupProlongation() override
	{}

	void SetupRestriction() override
	{
		double scalingFactor = 1.0;
		//scalingFactor = 1.0 / 4.0;
		R = scalingFactor * P.transpose();
	}

	void OnEndSetup() override
	{
		if (this->CoarserLevel)
		{
			CondensedLevel* coarse = dynamic_cast<CondensedLevel*>(this->CoarserLevel);
			coarse->A_T_T = &A_T_Tc;
			coarse->A_T_F = &A_T_Fc;
			coarse->A_F_F = &A_F_Fc;
			coarse->inv_A_T_T = inv_A_T_Tc;
		}
	}

	~CondensedLevel()
	{
		if (!this->IsFinestLevel())
			delete OperatorMatrix;
	}
};

class CondensedAMG : public Multigrid
{
private:
	CAMGProlongation _prolongation = CAMGProlongation::P;
	int _cellBlockSize;
	int _faceBlockSize;
public:

	CondensedAMG(int cellBlockSize, int faceBlockSize, CAMGProlongation prolongation, int nLevels = 0)
		: Multigrid(nLevels)
	{
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
		this->_prolongation = prolongation;
		this->BlockSizeForBlockSmoothers = faceBlockSize;
		this->UseGalerkinOperator = true;
		this->_fineLevel = new CondensedLevel(0, cellBlockSize, faceBlockSize, prolongation);
	}

	void BeginSerialize(ostream& os) const override
	{
		os << "CondensedAMG" << endl;
		os << "\t" << "Prolongation       : ";
		if (_prolongation == CAMGProlongation::P)
			os << "P ";
		else if (_prolongation == CAMGProlongation::P1P2)
			os << "P1P2 ";
		else if (_prolongation == CAMGProlongation::Q_F)
			os << "Q_F ";
		os << "[-prolong " << (unsigned)_prolongation << "]" << endl;
	}

	void EndSerialize(ostream& os) const override
	{
	}

	void Setup(const SparseMatrix& A) override
	{
		assert(false && "This Setup method cannot be used in this solver.");
	}

	void Setup(const SparseMatrix& A, const SparseMatrix& A_T_T, const SparseMatrix& A_T_F, const SparseMatrix& A_F_F) override
	{
		CondensedLevel* fine = dynamic_cast<CondensedLevel*>(this->_fineLevel);
		fine->A_T_T = &A_T_T;
		fine->A_T_F = &A_T_F;
		fine->A_F_F = &A_F_F;
		SparseMatrix* inv_A_T_T = new SparseMatrix(Utils::InvertBlockDiagMatrix(A_T_T, _cellBlockSize));
		fine->inv_A_T_T = inv_A_T_T;
		Multigrid::Setup(A);
	}

	Vector Solve(const Vector& b, string initialGuessCode) override
	{
		if (initialGuessCode.compare("smooth") == 0)
			Utils::Warning("Smooth initial guess unmanaged in AMG.");
		return Multigrid::Solve(b, initialGuessCode);
	}
	
protected:
	Level* CreateCoarseLevel(Level* fineLevel) override
	{
		CondensedLevel* coarseLevel = new CondensedLevel(fineLevel->Number + 1, _cellBlockSize, _faceBlockSize, _prolongation);
		return coarseLevel;
	}
};