#pragma once
#include <cstdio>
#include "Utils.h"
#include "ParallelLoop.h"
#include "../Mesh/Mesh.h"
#include "../FunctionalBasis/FunctionalBasis.h"
using namespace std;

class L2
{
public:
	template <int Dim>
	static double Error(Mesh<Dim>* mesh, const FunctionalBasis<Dim>& basis, const Vector& solution, DomFunction exactSolution)
	{
		cout << endl << "Computing L2 error..." << endl;
		
		struct ChunkResult
		{
			double absoluteError = 0;
			double normExactSolution = 0;
		};

		int degreeUsedForNonPolynomials = basis.GetDegree() + 2;

		ParallelLoop<Element<Dim>*, ChunkResult> parallelLoop(mesh->Elements);
		parallelLoop.Execute([&basis, &solution, exactSolution, &degreeUsedForNonPolynomials](Element<Dim>* element, ParallelChunk<ChunkResult>* chunk)
			{
				auto approximate = basis.GetApproximateFunction(solution, element->Number * basis.NumberOfLocalFunctionsInElement(element));
				chunk->Results.absoluteError += element->L2ErrorPow2(approximate, exactSolution, degreeUsedForNonPolynomials);
				chunk->Results.normExactSolution += element->Integral([exactSolution](DomPoint p) { return pow(exactSolution(p), 2); }, degreeUsedForNonPolynomials);
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
		return absoluteError / normExactSolution;
	}
};