#pragma once
#include <vector>
#include <map>
#include "../Utils/Utils.h"
#include "../Utils/DiffusionPartition.h"
#include "../FunctionalBasis/FunctionalBasis.h"

template <int Dim>
class Face;


template <int Dim>
class Element 
{
private:
	std::map<Face<Dim>*, int> _facesLocalNumbering;
public:
	BigNumber Number;
	std::vector<Face<Dim>*> Faces;

	std::vector<Element<Dim>*> FinerElements;
	Element<Dim>* CoarserElement;
	std::vector<Face<Dim>*> FinerFacesRemoved;

	Element(BigNumber number)
	{
		this->Number = number;
	}

	Element<Dim>* ElementOnTheOtherSideOf(Face<Dim>* face)
	{
		return face->GetNeighbour(this);
	}

	int LocalNumberOf(Element<Dim>* finerElement)
	{
		for (int i = 0; i < this->FinerElements.size(); i++)
		{
			if (this->FinerElements[i] == finerElement)
				return i;
		}
		assert(false);
	}

	int LocalNumberOf(Face<Dim>* face)
	{
		return this->_facesLocalNumbering[face];
	}

	virtual function<double(Point)> EvalPhiOnFace(Face<Dim>* face, BasisFunction<Dim>* phi)
	{
		function<double(Point)> evalOnFace = [this, face, phi](RefPoint refPoint1D) {
			DomPoint domainPoint2D = face->ConvertToDomain(refPoint1D);
			RefPoint refPoint2D = this->ConvertToReference(domainPoint2D);
			return phi->Eval(refPoint2D);
		};
		return evalOnFace;
	}

	virtual function<vector<double>(Point)> GradPhiOnFace(Face<Dim>* face, BasisFunction<Dim>* phi)
	{
		function<vector<double>(Point)> gradOnFace = [this, face, phi](RefPoint refPoint1D) {
			DomPoint domainPoint2D = face->ConvertToDomain(refPoint1D);
			RefPoint refPoint2D = this->ConvertToReference(domainPoint2D);
			return phi->Grad(refPoint2D);
		};
		return gradOnFace;
	}

	double FrontierMeasure()
	{
		double measure = 0;
		for (auto f : this->Faces)
			measure += f->Measure();
	}

	virtual double GetDiameter() = 0;
	virtual double Measure() = 0;
	virtual DomPoint ConvertToDomain(RefPoint refPoint) = 0;
	virtual RefPoint ConvertToReference(DomPoint domainPoint) = 0;
	virtual vector<double> OuterNormalVector(Face<Dim>* face) = 0;
	virtual double Integral(BasisFunction<Dim>* phi) = 0;
	virtual double IntegralGlobalFunction(function<double(DomPoint)> globalFunction) = 0;
	virtual double ComputeIntegral(function<double(RefPoint)> func, int numberOfDerivatives) = 0;
	virtual double ComputeIntegral(function<double(RefPoint)> func, int numberOfDerivatives, int polynomialDegree) = 0;
	virtual double L2ErrorPow2(function<double(RefPoint)> approximate, function<double(DomPoint)> exactSolution) = 0;
	virtual double DiffusionCoefficient(DiffusionPartition diffusionPartition) = 0;
	virtual vector<RefPoint> GetNodalPoints(FunctionalBasis<Dim>* basis) = 0;

	virtual ~Element() {}

protected:
	void AddFace(Face<Dim>* face)
	{
		this->Faces.push_back(face);

		int faceLocalNumber = static_cast<int>(this->_facesLocalNumbering.size());
		this->_facesLocalNumbering.insert(std::pair<Face<Dim>*, int>(face, faceLocalNumber));
	}
public:
	static double InnerProduct(vector<double> vector1, vector<double> vector2)
	{
		double innerProduct = 0;
		for (int i=0; i<Dim; i++)
			innerProduct += vector1[i] * vector2[i];
		return innerProduct;
	}
};