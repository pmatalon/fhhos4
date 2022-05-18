#pragma once
#include "../../Mesh/Mesh.h"
#include "../../Utils/Utils.h"
#include "../../TestCases/Diffusion/DiffusionTestCase.h"
#include "../../Geometry/CartesianShape.h"
#include "../../Geometry/2D/Triangle.h"
#include "../../Utils/ParallelLoop.h"
#include "../../Utils/ExportModule.h"
#include "Diff_DGElement.h"
#include "Diff_DGFace.h"
using namespace std;

template <int Dim>
class Diffusion_DG
{
public:
	Mesh<Dim>* _mesh;
	FunctionalBasis<Dim>* Basis;
	SparseMatrix A;
	Vector b;
	Vector SystemSolution;
private:
	DiffusionTestCase<Dim>* _testCase;
	bool _autoPenalization;
	int _penalizationCoefficient;
public:
	Diffusion_DG(Mesh<Dim>* mesh, DiffusionTestCase<Dim>* testCase, string outputDirectory, FunctionalBasis<Dim>* basis, int penalizationCoefficient)
	{ 
		this->_mesh = mesh;
		this->_testCase = testCase;
		this->Basis = basis;
		this->_autoPenalization = penalizationCoefficient == -1;
		if (_autoPenalization)
			this->_penalizationCoefficient = pow(Dim, 2) * pow(Basis->GetDegree() + 1, 2) / this->_mesh->H(); // Ralph-Hartmann
		else
			this->_penalizationCoefficient = penalizationCoefficient;

		//Problem<Dim>::AddFilePrefix("_DG_SIPG_" + basis->Name() + "_pen" + to_string(penalizationCoefficient));
	}

	void PrintDiscretization()
	{
		cout << "Mesh: " << this->_mesh->Description() << endl;
		cout << "Discretization: Discontinuous Galerkin SIPG" << endl;
		cout << "\tPolynomial space: " << (Basis->UsePolynomialSpaceQ() ? "Q" : "P") << endl;
		cout << "\tPolynomial basis: " << Basis->Name() << endl;
		cout << "\tPenalization coefficient: " << _penalizationCoefficient << (_autoPenalization ? " (automatic)" : "") << endl;
		cout << "Local functions: " << Basis->Size() << endl;
		for (BasisFunction<Dim>* phi : Basis->LocalFunctions())
			cout << "\t " << phi->ToString() << endl;
		BigNumber nUnknowns = static_cast<int>(this->_mesh->Elements.size()) * Basis->Size();
		cout << "Unknowns   : " << nUnknowns << endl;
	}

	double L2Error(DomFunction exactSolution)
	{
		struct ChunkResult
		{
			double absoluteError = 0;
			double normExactSolution = 0;
		};

		ParallelLoop<Element<Dim>*, ChunkResult> parallelLoop(_mesh->Elements);
		parallelLoop.Execute([this, exactSolution](Element<Dim>* element, ParallelChunk<ChunkResult>* chunk)
			{
				auto approximate = Basis->GetApproximateFunction(SystemSolution, element->Number * Basis->Size());
				chunk->Results.absoluteError += element->L2ErrorPow2(approximate, exactSolution);
				chunk->Results.normExactSolution += element->Integral([exactSolution](const DomPoint& p) { return pow(exactSolution(p), 2); });
			});


		double absoluteError = 0;
		double normExactSolution = 0;

		parallelLoop.AggregateChunkResults([&absoluteError, &normExactSolution](ChunkResult chunkResult)
			{
				absoluteError += chunkResult.absoluteError;
				normExactSolution += chunkResult.normExactSolution;
			});

		absoluteError = sqrt(absoluteError);
		normExactSolution = sqrt(normExactSolution);
		return normExactSolution != 0 ? absoluteError / normExactSolution : absoluteError;
	}

	void Assemble(const ActionsArguments& actions, const ExportModule& out) //override
	{
		auto mesh = this->_mesh;
		auto basis = this->Basis;
		auto penalizationCoefficient = this->_penalizationCoefficient;

		if (actions.LogAssembly)
			this->PrintDiscretization();


		BigNumber nUnknowns = static_cast<int>(mesh->Elements.size()) * basis->Size();
		this->b = Vector(nUnknowns);

		if (actions.LogAssembly)
		{
			cout << "--------------------------------------------------------" << endl;
			cout << "Assembly..." << endl;
		}

		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreMassMatrix(basis);
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreStiffnessMatrix(basis);
		if (Dim == 2)
			Triangle::InitReferenceShape()->ComputeAndStoreMassMatrix((FunctionalBasis<2>*)basis);
		
		//--------------------------------------------//
		// Iteration on the elements: diagonal blocks //
		//--------------------------------------------//

		ParallelLoop<Element<Dim>*, EmptyResultChunk> parallelLoop(mesh->Elements);

		vector<NonZeroCoefficients> chunksMatrixCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksMassMatrixCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksVolumicCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksCouplingCoeffs(parallelLoop.NThreads);
		vector<NonZeroCoefficients> chunksPenCoeffs(parallelLoop.NThreads);
		
		for (unsigned int threadNumber = 0; threadNumber < parallelLoop.NThreads; threadNumber++)
		{
			ParallelChunk<EmptyResultChunk>* chunk = parallelLoop.Chunks[threadNumber];

			chunk->ThreadFuture = std::async([this, mesh, basis, &actions, penalizationCoefficient, chunk, &chunksMatrixCoeffs, &chunksMassMatrixCoeffs, &chunksVolumicCoeffs, &chunksCouplingCoeffs, &chunksPenCoeffs]()
			{
				BigNumber nnzApproximate = chunk->Size() * basis->Size() * (2 * Dim + 1);
				NonZeroCoefficients matrixCoeffs(nnzApproximate);
				NonZeroCoefficients massMatrixCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);
				NonZeroCoefficients volumicCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);
				NonZeroCoefficients couplingCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);
				NonZeroCoefficients penCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);

				for (BigNumber iElem = chunk->Start; iElem < chunk->End; iElem++)
				{
					Diff_DGElement<Dim>* element = dynamic_cast<Diff_DGElement<Dim>*>(mesh->Elements[iElem]);
					//cout << "Element " << element->Number << endl;

					for (BasisFunction<Dim>* phi1 : basis->LocalFunctions())
					{
						BigNumber basisFunction1 = element->Number * basis->Size() + phi1->LocalNumber;

						// Current element (block diagonal)
						for (BasisFunction<Dim>* phi2 : basis->LocalFunctions())
						{
							BigNumber basisFunction2 = element->Number * basis->Size() + phi2->LocalNumber;

							//cout << "\t phi" << phi1->LocalNumber << " = " << phi1->ToString() << " phi" << phi2->LocalNumber << " = " << phi2->ToString() << endl;

							double volumicTerm = element->VolumicTerm(phi1, phi2);
							//cout << "\t\t volumic = " << volumicTerm << endl;

							double coupling = 0;
							double penalization = 0;
							for (Face<Dim>* f : element->Faces)
							{
								Diff_DGFace<Dim>* face = dynamic_cast<Diff_DGFace<Dim>*>(f);

								double c = face->CouplingTerm(element, phi1, element, phi2);
								double p = face->PenalizationTerm(element, phi1, element, phi2, penalizationCoefficient);
								coupling += c;
								penalization += p;
								//cout << "\t\t " << face->ToString() << ":\t c=" << c << "\tp=" << p << endl;
							}

							//cout << "\t\t TOTAL = " << volumicTerm + coupling + penalization << endl;

							if (actions.Export.AssemblyTermMatrices)
							{
								volumicCoeffs.Add(basisFunction1, basisFunction2, volumicTerm);
								couplingCoeffs.Add(basisFunction1, basisFunction2, coupling);
								penCoeffs.Add(basisFunction1, basisFunction2, penalization);
							}
							matrixCoeffs.Add(basisFunction1, basisFunction2, volumicTerm + coupling + penalization);
							if (actions.Export.AssemblyTermMatrices)
							{
								double massTerm = element->MassTerm(phi1, phi2);
								massMatrixCoeffs.Add(basisFunction1, basisFunction2, massTerm);
							}
						}

						double rhs = element->SourceTerm(phi1, _testCase->SourceFunction);
						this->b(basisFunction1) = rhs;
					}
				}

				chunksMatrixCoeffs[chunk->ThreadNumber] = matrixCoeffs;
				chunksMassMatrixCoeffs[chunk->ThreadNumber] = massMatrixCoeffs;
				chunksVolumicCoeffs[chunk->ThreadNumber] = volumicCoeffs;
				chunksCouplingCoeffs[chunk->ThreadNumber] = couplingCoeffs;
				chunksPenCoeffs[chunk->ThreadNumber] = penCoeffs;
			}
			);
		}

		BigNumber nnzApproximate = mesh->Elements.size() * basis->Size() * (2 * Dim + 1);
		NonZeroCoefficients matrixCoeffs(nnzApproximate);
		NonZeroCoefficients massMatrixCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);
		NonZeroCoefficients volumicCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);
		NonZeroCoefficients couplingCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);
		NonZeroCoefficients penCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);

		parallelLoop.Wait();

		for (unsigned int threadNumber = 0; threadNumber < parallelLoop.NThreads; threadNumber++)
		{
			matrixCoeffs.Add(chunksMatrixCoeffs[threadNumber]);
			massMatrixCoeffs.Add(chunksMassMatrixCoeffs[threadNumber]);
			volumicCoeffs.Add(chunksVolumicCoeffs[threadNumber]);
			couplingCoeffs.Add(chunksCouplingCoeffs[threadNumber]);
			penCoeffs.Add(chunksPenCoeffs[threadNumber]);
		}

		chunksMatrixCoeffs.clear();
		chunksMassMatrixCoeffs.clear();
		chunksVolumicCoeffs.clear();
		chunksCouplingCoeffs.clear();
		chunksPenCoeffs.clear();

		//---------------------------------------------//
		// Iteration on the faces: off-diagonal blocks //
		//---------------------------------------------//

		ParallelLoop<Face<Dim>*, EmptyResultChunk> parallelLoopFaces(mesh->Faces);

		chunksMatrixCoeffs = vector<NonZeroCoefficients>(parallelLoopFaces.NThreads);
		chunksCouplingCoeffs = vector<NonZeroCoefficients>(parallelLoopFaces.NThreads);
		chunksPenCoeffs = vector<NonZeroCoefficients>(parallelLoopFaces.NThreads);


		for (unsigned int threadNumber = 0; threadNumber < parallelLoopFaces.NThreads; threadNumber++)
		{
			ParallelChunk<EmptyResultChunk>* chunk = parallelLoopFaces.Chunks[threadNumber];

			chunk->ThreadFuture = std::async([this, mesh, basis, &actions, penalizationCoefficient, chunk, &chunksMatrixCoeffs, &chunksCouplingCoeffs, &chunksPenCoeffs]()
			{
				BigNumber nnzApproximate = chunk->Size() * basis->Size() * (2 * Dim + 1);
				NonZeroCoefficients matrixCoeffs(nnzApproximate);
				NonZeroCoefficients couplingCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);
				NonZeroCoefficients penCoeffs(actions.Export.AssemblyTermMatrices ? nnzApproximate : 0);

				for (BigNumber iElem = chunk->Start; iElem < chunk->End; ++iElem)
				{
					Diff_DGFace<Dim>* face = dynamic_cast<Diff_DGFace<Dim>*>(mesh->Faces[iElem]);
					if (face->IsDomainBoundary)
						continue;

					//cout << "Face " << face->Number << endl;

					for (BasisFunction<Dim>* phi1 : basis->LocalFunctions())
					{
						BigNumber basisFunction1 = face->Element1->Number * basis->Size() + phi1->LocalNumber;
						for (BasisFunction<Dim>* phi2 : basis->LocalFunctions())
						{
							//cout << "\t phi" << phi1->LocalNumber << " = " << phi1->ToString() << " phi" << phi2->LocalNumber << " = " << phi2->ToString() << endl;

							BigNumber basisFunction2 = face->Element2->Number * basis->Size() + phi2->LocalNumber;
							double coupling = face->CouplingTerm(face->Element1, phi1, face->Element2, phi2);
							double penalization = face->PenalizationTerm(face->Element1, phi1, face->Element2, phi2, penalizationCoefficient);

							//cout << "\t\t\t c=" << coupling << "\tp=" << penalization << endl;

							if (actions.Export.AssemblyTermMatrices)
							{
								couplingCoeffs.Add(basisFunction1, basisFunction2, coupling);
								couplingCoeffs.Add(basisFunction2, basisFunction1, coupling);

								penCoeffs.Add(basisFunction1, basisFunction2, penalization);
								penCoeffs.Add(basisFunction2, basisFunction1, penalization);
							}
							matrixCoeffs.Add(basisFunction1, basisFunction2, coupling + penalization);
							matrixCoeffs.Add(basisFunction2, basisFunction1, coupling + penalization);
						}
					}
				}

				chunksMatrixCoeffs[chunk->ThreadNumber] = matrixCoeffs;
				chunksCouplingCoeffs[chunk->ThreadNumber] = couplingCoeffs;
				chunksPenCoeffs[chunk->ThreadNumber] = penCoeffs;
			}
			);
		}

		parallelLoopFaces.Wait();

		for (unsigned int threadNumber = 0; threadNumber < parallelLoopFaces.NThreads; threadNumber++)
		{
			matrixCoeffs.Add(chunksMatrixCoeffs[threadNumber]);
			couplingCoeffs.Add(chunksCouplingCoeffs[threadNumber]);
			penCoeffs.Add(chunksPenCoeffs[threadNumber]);
		}

		//---------------//
		// Matrix export //
		//---------------//

		this->A = SparseMatrix(nUnknowns, nUnknowns);
		matrixCoeffs.Fill(this->A);
		cout << "nnz(A) = " << this->A.nonZeros() << endl;

		if (actions.Export.LinearSystem)
		{
			cout << "Export of the linear system..." << endl;
			out.ExportMatrix(this->A, "A");

			out.ExportVector(this->b, "b");
		}

		if (actions.Export.AssemblyTermMatrices)
		{
			SparseMatrix M(nUnknowns, nUnknowns);
			massMatrixCoeffs.Fill(M);
			out.ExportMatrix(M, "Mass");

			SparseMatrix V(nUnknowns, nUnknowns);
			volumicCoeffs.Fill(V);
			out.ExportMatrix(V, "A_volumic");

			SparseMatrix C(nUnknowns, nUnknowns);
			couplingCoeffs.Fill(C);
			out.ExportMatrix(C, "A_coupling");

			SparseMatrix P(nUnknowns, nUnknowns);
			penCoeffs.Fill(P);
			out.ExportMatrix(P, "A_pen");
		}

	}
};

