#pragma once
#include "../../DG/Poisson_DG_Face.h"
#include "../../HHO/Poisson_HHO_Face.h"
#include "../CartesianShape.h"

class Edge : public Poisson_DG_Face<2>, public Poisson_HHO_Face<2>
{
private:
	double _width;
	DomPoint _center;

public:
	Vertex* Vertex1;
	Vertex* Vertex2;

	Edge(BigNumber number, Vertex* v1, Vertex* v2, Element<2>* element1, Element<2>* element2) : 
		Poisson_DG_Face(number, element1, element2),
		Poisson_HHO_Face(number, element1, element2),
		Face(number, element1, element2)
	{
		Vertex1 = v1;
		Vertex2 = v2;
		Init();
	}

	Edge(BigNumber number, Vertex* v1, Vertex* v2, Element<2>* element1) :
		Edge(number, v1, v2, element1, nullptr)
	{ }

	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	void Serialize(ostream& os) const override
	{
		Face<2>::Serialize(os);
		os << ": ";
		Vertex1->Serialize(os, 2);
		os << "--";
		Vertex2->Serialize(os, 2);
	}

	inline void Init()
	{
		_width = sqrt(pow(Vertex2->X - Vertex1->X, 2) + pow(Vertex2->Y - Vertex1->Y, 2));
		_center = DomPoint((Vertex1->X + Vertex2->X) / 2, (Vertex1->Y + Vertex2->Y) / 2);
	}

	inline double GetDiameter()
	{
		return _width;
	}
	inline double Measure()
	{
		return _width;
	}
	inline DomPoint Center()
	{
		return _center;
	}

	DomPoint ConvertToDomain(RefPoint referenceElementPoint)
	{
		double x1 = Vertex1->X;
		double x2 = Vertex2->X;
		double y1 = Vertex1->Y;
		double y2 = Vertex2->Y;

		double t = referenceElementPoint.X;

		DomPoint p((x2 - x1) / 2 * t + (x2 + x1) / 2, (y2 - y1) / 2 * t + (y2 + y1) / 2);
		return p;
	}

	RefPoint ConvertToReference(DomPoint domainPoint)
	{
		double x1 = Vertex1->X;
		double x2 = Vertex2->X;

		double x = domainPoint.X;

		RefPoint p(2 / (x2 - x1) * x - (x2 + x1) / (x2 - x1));
		return p;
	}

	double Integral(RefFunction func) const
	{
		double integralOnReferenceShape = ReferenceEdge().Integral(func);
		return DetJacobian() * integralOnReferenceShape;
	}

	double Integral(RefFunction func, int polynomialDegree) const
	{
		double integralOnReferenceShape = ReferenceEdge().Integral(func, polynomialDegree);
		return DetJacobian() * integralOnReferenceShape;
	}

	Face<2>* CreateSameGeometricFace(BigNumber number, Element<2>* element1)
	{
		Face<2>* copy = new Edge(number, Vertex1, Vertex2, element1);
		copy->IsDomainBoundary = this->IsDomainBoundary;
		return copy;
	}

	void ExportFaceToMatlab(FILE* file)
	{
		//fprintf(file, "%llu %.17g %.17g %.17g %.17g %.17g %.17g %d %d\n", this->Number, this->Origin->X, this->Origin->Y, this->Origin->Z, this->WidthX, this->WidthY, this->WidthZ, this->Orientation, this->IsDomainBoundary);
	}

private:
	inline static ReferenceCartesianShape<1> ReferenceEdge()
	{
		return CartesianShape<2, 1>::ReferenceShape;
	}

	inline double DetJacobian() const
	{
		return _width / 2;
	}

	//----------------------------------------------------------------//
	//                 Poisson_HHO_Face implementation                //
	//----------------------------------------------------------------//

public:
	DenseMatrix FaceMassMatrix(FunctionalBasis<1>* basis)
	{
		return DetJacobian() * CartesianShape<2, 1>::ReferenceShape.FaceMassMatrix(basis);
	}
};