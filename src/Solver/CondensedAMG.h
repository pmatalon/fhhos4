#pragma once
#include <mutex>
#include "HybridAlgebraicMesh.h"
#include "AlgebraicMesh.h"
#include "Multigrid.h"
using namespace std;

class CondensedLevel : public Level
{
private:
	CAMGFaceProlongation _faceProlong = CAMGFaceProlongation::FaceAggregates;
	CAMGProlongation _multigridProlong = CAMGProlongation::ReconstructionTrace1Step;
	int _cellBlockSize;
	int _faceBlockSize;
	double _strongCouplingThreshold;
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
	CondensedLevel(int number, int cellBlockSize, int faceBlockSize, double strongCouplingThreshold, CAMGFaceProlongation faceProlong, CAMGProlongation prolongation)
		: Level(number)
	{
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
		this->_strongCouplingThreshold = strongCouplingThreshold;
		this->_multigridProlong = prolongation;
		this->_faceProlong = faceProlong;
	}

	BigNumber NUnknowns() override
	{
		return A_F_F->rows();
	}

	void ExportVector(const Vector& v, string suffix, int levelNumber) override
	{ }

	void ExportMatrix(const SparseMatrix& M, string suffix, int levelNumber) override
	{ 
		Eigen::saveMarket(M, Utils::ProgramArgs.OutputDirectory + "/" + suffix + ".dat");
	}

	void CoarsenMesh(CoarseningStrategy coarseningStgy, int coarseningFactor, bool& noCoarserMeshProvided, bool& coarsestPossibleMeshReached) override
	{
		noCoarserMeshProvided = false;
		coarsestPossibleMeshReached = false;

		if (!this->OperatorMatrix)
		{
			if (this->UseGalerkinOperator)
				ComputeGalerkinOperator();
			else
				SetupDiscretizedOperator();
		}

		const SparseMatrix* A_T_T1 = A_T_T;
		const SparseMatrix* A_T_F1 = A_T_F;
		const SparseMatrix* A_F_F1 = A_F_F;
		SparseMatrix *A_T_T2, *A_T_F2, *A_F_F2, *P, *Q_F;
		double coarseningRatio = 0;
		double nCoarsenings = 0;
		//while (coarseningRatio < 3)
		for (int i = 0; i < 2; i++)
		{
			cout << "\tPairwise aggregation" << endl;

			std::tie(A_T_T2, A_T_F2, A_F_F2, P, Q_F, coarsestPossibleMeshReached) = PairwiseAggregate(*A_T_T1, *A_T_F1, *A_F_F1, coarseningStgy);
			if (coarsestPossibleMeshReached)
				return;

			if (nCoarsenings > 0)
				delete A_T_T1, A_T_F1, A_F_F1;

			A_T_T1 = A_T_T2;
			A_T_F1 = A_T_F2;
			A_F_F1 = A_F_F2;

			if (nCoarsenings == 0)
			{
				this->P = *P;
				this->Q_F = *Q_F;
			}
			else
			{
				this->P = this->P * *P;
				this->Q_F = this->Q_F * *Q_F;
			}

			delete P, Q_F;

			double nFine = this->A_T_F->cols();
			double nCoarse = A_T_F1->cols();
			coarseningRatio = nFine / nCoarse;
			nCoarsenings++;
		}

		this->A_T_Tc = *A_T_T2;
		this->A_T_Fc = *A_T_F2;
		this->A_F_Fc = *A_F_F2;

		/*
		// Assembly of the double-aggregates
		vector<vector<HybridAlgebraicElement*>> doubleElemAggregates(coarseMesh._coarseElements.size());
		NumberParallelLoop<EmptyResultChunk> parallelLoopE(coarseMesh._coarseElements.size());
		parallelLoopE.Execute([this, &coarseMesh, &doubleElemAggregates](BigNumber secondAggregNumber)
			{
				HybridElementAggregate& agg2 = coarseMesh._coarseElements[secondAggregNumber];
				vector<HybridAlgebraicElement*> doubleAggregate;
				for (HybridAlgebraicElement* ce1 : agg2.FineElements)
				{
					HybridElementAggregate& agg1 = mesh._coarseElements[ce1->Number];
					for (HybridAlgebraicElement* fe : agg1.FineElements)
						doubleAggregate.push_back(fe);
				}
				doubleElemAggregates[secondAggregNumber] = doubleAggregate;
			});

		vector<vector<HybridAlgebraicFace*>> doubleFaceAggregates(coarseMesh._coarseFaces.size());
		NumberParallelLoop<EmptyResultChunk> parallelLoopF(coarseMesh._coarseFaces.size());
		parallelLoopF.Execute([this, &coarseMesh, &doubleFaceAggregates](BigNumber secondAggregNumber)
			{
				HybridFaceAggregate& agg2 = coarseMesh._coarseFaces[secondAggregNumber];
				vector<HybridAlgebraicFace*> doubleAggregate;
				for (HybridAlgebraicFace* cf1 : agg2.FineFaces)
				{
					HybridFaceAggregate& agg1 = mesh._coarseFaces[cf1->Number];
					for (HybridAlgebraicFace* fe : agg1.FineFaces)
						doubleAggregate.push_back(fe);
				}
				doubleFaceAggregates[secondAggregNumber] = doubleAggregate;
			});
		
		//SparseMatrix Q_T = BuildQ_T(doubleElemAggregates);
		//SparseMatrix FaceProlongation = BuildQ_F(doubleFaceAggregates);
		*/

	}

	// Returns <A_T_Tc, A_T_Fc, A_F_Fc, P, Q_F, coarsestPossibleMeshReached>
	tuple<SparseMatrix*, SparseMatrix*, SparseMatrix*, SparseMatrix*, SparseMatrix*, bool> PairwiseAggregate(const SparseMatrix A_T_T, const SparseMatrix A_T_F, const SparseMatrix A_F_F, CoarseningStrategy coarseningStgy)
	{
		HybridAlgebraicMesh mesh(_cellBlockSize, _faceBlockSize, _strongCouplingThreshold);
		mesh.Build(A_T_T, A_T_F, A_F_F);

		bool coarsestPossibleMeshReached = false;
		mesh.PairWiseAggregate(coarseningStgy, coarsestPossibleMeshReached);
		if (coarsestPossibleMeshReached)
			return { nullptr, nullptr, nullptr, nullptr, nullptr, coarsestPossibleMeshReached };

		/*cout << "------------- Elem" << endl;
		for (HybridElementAggregate& agg : mesh._coarseElements)
		{
			cout << "(";
			for (auto e : agg.FineElements)
				cout << e->Number << ", ";
			cout << ")" << endl;
		}

		cout << "------------- Faces" << endl;
		for (HybridFaceAggregate& agg : mesh._coarseFaces)
		{
			cout << "(";
			for (auto f : agg.FineFaces)
				cout << f->Number << ", ";
			cout << ")" << endl;
		}*/

		// Cell-prolongation operator with only one 1 coefficient per row
		SparseMatrix Q_T = BuildQ_T(mesh);

		SparseMatrix* Q_F;
		if (this->_faceProlong == CAMGFaceProlongation::BoundaryAggregatesInteriorAverage)
		{
			// Face-prolongation operator with only one 1 coefficient per row for kept or aggregated faces, average for removed faces
			Q_F = new SparseMatrix(BuildQ_F(mesh, A_F_F));
		}
		else if (this->_faceProlong == CAMGFaceProlongation::BoundaryAggregatesInteriorZero)
		{
			Q_F = new SparseMatrix(BuildQ_F_0Interior(mesh));
		}
		else if (this->_faceProlong == CAMGFaceProlongation::FaceAggregates)
		{
			if (coarseningStgy != CoarseningStrategy::CAMGAggregFaces)
				Utils::FatalError("This coarsening strategy is incompatible with this face prolongation operator.");

			//ExportMatrix(*A_F_F, "A_F_F", 0);
			AlgebraicMesh skeleton(_faceBlockSize, 0);
			//skeleton.Build(*A_F_F);
			skeleton.Build(*this->OperatorMatrix);
			skeleton.PairWiseAggregate(coarsestPossibleMeshReached);
			if (coarsestPossibleMeshReached)
				return { nullptr, nullptr, nullptr, nullptr, nullptr, coarsestPossibleMeshReached };
			Q_F = new SparseMatrix(BuildQ_F_AllAggregated(skeleton));
		}
		else
			Utils::FatalError("Unmanaged -face-prolong");

		// Intermediate coarse operators
		SparseMatrix* A_T_Tc     = new SparseMatrix(Q_T.transpose() * A_T_T * Q_T);
		SparseMatrix* A_T_Fc_tmp = new SparseMatrix(Q_T.transpose() * A_T_F * *Q_F);

		/*ExportMatrix(A_T_F, "A_T_F", 0);
		ExportMatrix(Q_T, "Q_T", 0);
		ExportMatrix(*Q_F, "Q_F", 0);
		ExportMatrix(*A_T_Tc, "A_T_Tc", 0);
		ExportMatrix(*A_T_Fc_tmp, "A_T_Fc", 0);*/

		SparseMatrix* P;
		if (_multigridProlong == CAMGProlongation::ReconstructionTrace1Step || _multigridProlong == CAMGProlongation::ReconstructionTrace2Steps) // 1
		{
			SparseMatrix inv_A_T_Tc = Utils::InvertBlockDiagMatrix(*A_T_Tc, _cellBlockSize);
			// Theta: reconstruction from the coarse faces to the coarse cells
			SparseMatrix Theta = -inv_A_T_Tc * *A_T_Fc_tmp;
			// Pi: average on both sides of each face
			SparseMatrix Pi = BuildTrace(mesh);

			P = new SparseMatrix(Pi * Q_T * Theta);
		}
		else if (_multigridProlong == CAMGProlongation::FaceProlongation) // 3
		{
			// -g 1 -prolong 3 -face-prolong 3 -cs z
			P = Q_F;
		}
		else if (_multigridProlong == CAMGProlongation::ReconstructionTranspose2Steps) // 4
		{
			/*SparseMatrix inv_A_T_T1 = Utils::InvertBlockDiagMatrix(A_T_T1, _cellBlockSize);
			// Theta: reconstruction from the coarse faces to the coarse cells
			SparseMatrix Theta1 = -inv_A_T_T1 * A_T_F1;
			// Pi: transpose of Theta
			SparseMatrix inv_A_T_T = Utils::InvertBlockDiagMatrix(*A_T_T, _cellBlockSize);
			SparseMatrix Theta = -inv_A_T_T * *A_T_F;
			SparseMatrix P1 = Theta.transpose() * Q_T1 * Theta1;

			this->inv_A_T_Tc = new SparseMatrix(Utils::InvertBlockDiagMatrix(A_T_Tc, _cellBlockSize));
			SparseMatrix Theta2 = -(*inv_A_T_Tc) * A_T_Fc;
			SparseMatrix P2 = Theta1.transpose() * Q_T2 * Theta2;
			this->P = P1 * P2;*/
		}
		else if (_multigridProlong == CAMGProlongation::ReconstructTraceOrInject) // 5
		{
			SparseMatrix inv_A_T_Tc = Utils::InvertBlockDiagMatrix(*A_T_Tc, _cellBlockSize);
			SparseMatrix Theta = -inv_A_T_Tc * *A_T_Fc_tmp;   // Reconstruct
			SparseMatrix Pi = BuildTraceOnRemovedFaces(mesh); // Trace
												              // Inject = Q_F
			SparseMatrix ReconstructAndTrace1 = Pi * Q_T * Theta;

			NumberParallelLoop<CoeffsChunk> parallelLoop(mesh._faces.size());
			parallelLoop.ReserveChunkCoeffsSize(_faceBlockSize * 2 * _faceBlockSize);
			parallelLoop.Execute([this, &ReconstructAndTrace1, &Q_F, &mesh](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
				{
					const HybridAlgebraicFace* face = &mesh._faces[faceNumber];
					if (face->IsRemovedOnCoarseMesh)
						chunk->Results.Coeffs.CopyRows(faceNumber*_faceBlockSize, _faceBlockSize, ReconstructAndTrace1);
					else
						chunk->Results.Coeffs.CopyRows(faceNumber*_faceBlockSize, _faceBlockSize, *Q_F);
				});
			P = new SparseMatrix(Q_F->rows(), Q_F->cols());
			parallelLoop.Fill(*P);
		}
		else
			Utils::FatalError("Unmanaged prolongation");

		//SparseMatrix* A_T_Fc = new SparseMatrix(Q_T.transpose() * A_T_F * *P); // Kills -prolong 1 or 2
		SparseMatrix* A_T_Fc = A_T_Fc_tmp;
		SparseMatrix* A_F_Fc = new SparseMatrix(P->transpose() * A_F_F * *P);

		return { A_T_Tc, A_T_Fc, A_F_Fc, P, Q_F, coarsestPossibleMeshReached };
	}

	void SetupDiscretizedOperator() override
	{
		//SparseMatrix* schur = new SparseMatrix(*A_F_F - (A_T_F->transpose()) * (*inv_A_T_T) * (*A_T_F));
		//this->OperatorMatrix = schur;

		CondensedLevel* fine = dynamic_cast<CondensedLevel*>(this->FinerLevel);
		this->OperatorMatrix = new SparseMatrix(fine->Q_F.transpose() * *(fine->OperatorMatrix) * fine->Q_F);
	}

private:
	// Cell prolongation Q_T with only one 1 coefficient per row
	SparseMatrix BuildQ_T(const HybridAlgebraicMesh& mesh)
	{
		DenseMatrix Id = DenseMatrix::Identity(_cellBlockSize, _cellBlockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoopQ_T(mesh._elements.size());
		parallelLoopQ_T.Execute([this, &mesh, &Id](BigNumber elemNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridAlgebraicElement& elem = mesh._elements[elemNumber];
				chunk->Results.Coeffs.Add(elem.Number*_cellBlockSize, elem.CoarseElement->Number*_cellBlockSize, Id);
			});
		SparseMatrix Q_T = SparseMatrix(mesh._elements.size()*_cellBlockSize, mesh._coarseElements.size()*_cellBlockSize);
		parallelLoopQ_T.Fill(Q_T);
		return Q_T;
	}

	// Face prolongation FaceProlongation
	SparseMatrix BuildQ_F(const HybridAlgebraicMesh& mesh, const SparseMatrix A_F_F)
	{
		bool enableAnisotropyManagement = false;
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);

		NumberParallelLoop<CoeffsChunk> parallelLoop1(mesh._faces.size());
		parallelLoop1.Execute([this, &mesh, &Id, enableAnisotropyManagement, &A_F_F](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridAlgebraicFace* face = &mesh._faces[faceNumber];
				if (face->IsRemovedOnCoarseMesh && (!Utils::ProgramArgs.Solver.MG.ManageAnisotropy || !face->CoarseFace))
				{
					// Take the average value of the coarse element faces
					HybridElementAggregate* elemAggreg = face->Elements[0]->CoarseElement;

					if (enableAnisotropyManagement)
					{
						map<HybridFaceAggregate*, double> couplings;
						double totalCouplings = 0;
						for (HybridFaceAggregate* coarseFace : elemAggreg->CoarseFaces)
						{
							double avgCouplingCoarseFace = 0;
							for (HybridAlgebraicFace* f : coarseFace->FineFaces)
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
							HybridFaceAggregate* coarseFace = it->first;
							double coupling = it->second;
							chunk->Results.Coeffs.Add(face->Number, coarseFace->Number, -coupling / abs(totalCouplings) * Id);
						}
					}
					else
					{
						for (HybridFaceAggregate* coarseFace : elemAggreg->CoarseFaces)
							chunk->Results.Coeffs.Add(face->Number*_faceBlockSize, coarseFace->Number*_faceBlockSize, 1.0 / elemAggreg->CoarseFaces.size() * Id);
					}
				}
				else
					chunk->Results.Coeffs.Add(face->Number*_faceBlockSize, face->CoarseFace->Number*_faceBlockSize, Id);
			});


		SparseMatrix Q_F = SparseMatrix(mesh._faces.size()*_faceBlockSize, mesh._coarseFaces.size()*_faceBlockSize);
		parallelLoop1.Fill(Q_F);
		return Q_F;
	}

	SparseMatrix BuildQ_F_0Interior(const HybridAlgebraicMesh& mesh)
	{
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);

		NumberParallelLoop<CoeffsChunk> parallelLoop(mesh._faces.size());
		parallelLoop.Execute([this, &mesh, &Id](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridAlgebraicFace* face = &mesh._faces[faceNumber];
				if (!face->IsRemovedOnCoarseMesh)
					chunk->Results.Coeffs.Add(face->Number*_faceBlockSize, face->CoarseFace->Number*_faceBlockSize, Id);
			});


		SparseMatrix Q_F = SparseMatrix(mesh._faces.size()*_faceBlockSize, mesh._coarseFaces.size()*_faceBlockSize);
		parallelLoop.Fill(Q_F);
		return Q_F;
	}

	SparseMatrix BuildQ_F_AllAggregated(const AlgebraicMesh& skeleton)
	{
		DenseMatrix Id = DenseMatrix::Identity(_faceBlockSize, _faceBlockSize);
		NumberParallelLoop<CoeffsChunk> parallelLoop(skeleton._elements.size());
		parallelLoop.Execute([this, &skeleton, &Id](BigNumber elemNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const AlgebraicElement& elem = skeleton._elements[elemNumber];
				chunk->Results.Coeffs.Add(elem.Number*_faceBlockSize, elem.CoarseElement->Number*_faceBlockSize, Id);
			});
		SparseMatrix Q_F = SparseMatrix(skeleton._elements.size()*_faceBlockSize, skeleton._coarseElements.size()*_faceBlockSize);
		parallelLoop.Fill(Q_F);
		return Q_F;
	}

	SparseMatrix BuildTrace(const HybridAlgebraicMesh& mesh)
	{
		// Pi: average on both sides of each face
		DenseMatrix traceOfConstant = DenseMatrix::Zero(_faceBlockSize, _cellBlockSize);
		traceOfConstant(0, 0) = 1;

		NumberParallelLoop<CoeffsChunk> parallelLoopPi(mesh._faces.size());
		parallelLoopPi.Execute([this, &mesh, &traceOfConstant](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridAlgebraicFace& face = mesh._faces[faceNumber];
				if (face.IsRemovedOnCoarseMesh)
					chunk->Results.Coeffs.Add(faceNumber*_faceBlockSize, face.Elements[0]->Number*_cellBlockSize, traceOfConstant);
				else
				{
					assert(!face.Elements.empty());
					for (HybridAlgebraicElement* elem : face.Elements)
						chunk->Results.Coeffs.Add(faceNumber*_faceBlockSize, elem->Number*_cellBlockSize, 1.0 / face.Elements.size()*traceOfConstant);
				}
			});
		SparseMatrix Pi = SparseMatrix(mesh._faces.size()*_faceBlockSize, mesh._elements.size()*_cellBlockSize);
		parallelLoopPi.Fill(Pi);
		return Pi;
	}

	SparseMatrix BuildTraceOnRemovedFaces(const HybridAlgebraicMesh& mesh)
	{
		// Pi: average on both sides of each face
		DenseMatrix traceOfConstant = DenseMatrix::Zero(_faceBlockSize, _cellBlockSize);
		traceOfConstant(0, 0) = 1;

		NumberParallelLoop<CoeffsChunk> parallelLoopPi(mesh._faces.size());
		parallelLoopPi.Execute([this, &mesh, &traceOfConstant](BigNumber faceNumber, ParallelChunk<CoeffsChunk>* chunk)
			{
				const HybridAlgebraicFace& face = mesh._faces[faceNumber];
				if (face.IsRemovedOnCoarseMesh)
					chunk->Results.Coeffs.Add(faceNumber*_faceBlockSize, face.Elements[0]->Number*_cellBlockSize, traceOfConstant);
			});
		SparseMatrix Pi = SparseMatrix(mesh._faces.size()*_faceBlockSize, mesh._elements.size()*_cellBlockSize);
		parallelLoopPi.Fill(Pi);
		return Pi;
	}

public:
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
		if (!this->IsFinestLevel() && !this->UseGalerkinOperator)
			delete OperatorMatrix;
	}
};

class CondensedAMG : public Multigrid
{
private:
	CAMGFaceProlongation _faceProlong = CAMGFaceProlongation::FaceAggregates;
	CAMGProlongation _multigridProlong = CAMGProlongation::ReconstructionTrace1Step;
	int _cellBlockSize;
	int _faceBlockSize;
	double _strongCouplingThreshold;
public:

	CondensedAMG(int cellBlockSize, int faceBlockSize, double strongCouplingThreshold, CAMGFaceProlongation faceProlong, CAMGProlongation prolongation, int nLevels = 0)
		: Multigrid(nLevels)
	{
		this->_cellBlockSize = cellBlockSize;
		this->_faceBlockSize = faceBlockSize;
		this->_strongCouplingThreshold = strongCouplingThreshold;
		this->_multigridProlong = prolongation;
		this->_faceProlong = faceProlong;
		this->BlockSizeForBlockSmoothers = faceBlockSize;
		this->UseGalerkinOperator = true;
		this->_fineLevel = new CondensedLevel(0, cellBlockSize, faceBlockSize, strongCouplingThreshold, faceProlong, prolongation);
	}

	void BeginSerialize(ostream& os) const override
	{
		os << "CondensedAMG" << endl;
		os << "\t" << "Prolongation       : ";
		if (_multigridProlong == CAMGProlongation::ReconstructionTrace1Step)
			os << "ReconstructionTrace1Step ";
		else if (_multigridProlong == CAMGProlongation::ReconstructionTrace2Steps)
			os << "ReconstructionTrace2Steps ";
		else if (_multigridProlong == CAMGProlongation::FaceProlongation)
			os << "FaceProlongation ";
		os << "[-prolong " << (unsigned)_multigridProlong << "]" << endl;
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
		CondensedLevel* coarseLevel = new CondensedLevel(fineLevel->Number + 1, _cellBlockSize, _faceBlockSize, _strongCouplingThreshold, _faceProlong, _multigridProlong);
		return coarseLevel;
	}
};