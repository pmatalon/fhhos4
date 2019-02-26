#pragma once
#include <vector>
#include "../Utils/Utils.h"
#include "../Utils/DiffusionPartition.h"
#include "../FunctionalBasis/BasisFunction.h"

template <short Dim>
class Face;

enum class StandardElementCode
{
	None,
	Interval,
	Square,
	Cube
};

template <short Dim>
class Element 
{
public:
	BigNumber Number;
	std::vector<Face<Dim>*> Faces;

	Element(BigNumber number)
	{
		this->Number = number;
	}

	virtual StandardElementCode StdElementCode() = 0;

	virtual double* OuterNormalVector(Face<Dim>* face) = 0;

	virtual double Integral(function<double(Point)> func) = 0;
	virtual function<double(Point)> EvalPhiOnFace(Face<Dim>* face, BasisFunction<Dim>* phi) = 0;
	virtual function<double*(Point)> GradPhiOnFace(Face<Dim>* face, BasisFunction<Dim>* phi) = 0;
	virtual double L2ErrorPow2(function<double(Point)> approximate, function<double(Point)> exactSolution) = 0;
	virtual double DiffusionCoefficient(DiffusionPartition diffusionPartition) = 0;

	virtual ~Element() {}
};