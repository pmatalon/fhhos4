#pragma once
#include "Vertex.h"
#include "Face.h"
#include "CartesianShape.h"
#include "../HHO/Poisson_HHO_Face.h"

template <int Dim>
class CartesianFace : public Poisson_DG_Face<Dim>, public Poisson_HHO_Face<Dim>, public CartesianShape<Dim, Dim-1>
{
public:

	CartesianFace(BigNumber number, Vertex* origin, double width, Element<Dim>* element1, Element<Dim>* element2, CartesianShapeOrientation orientation) :
		Poisson_DG_Face<Dim>(number, element1, element2),
		Poisson_HHO_Face<Dim>(number, element1, element2), 
		CartesianShape<Dim, Dim - 1>(origin, width, orientation)
	{}

	CartesianFace(BigNumber number, Vertex* origin, double firstWidth, double secondWidth, Element<Dim>* element1, Element<Dim>* element2, CartesianShapeOrientation orientation) :
		Poisson_DG_Face<Dim>(number, element1, element2),
		Poisson_HHO_Face<Dim>(number, element1, element2),
		CartesianShape<Dim, Dim - 1>(origin, firstWidth, secondWidth, orientation)
	{
		assert(Dim == 3);
	}

private:
	CartesianFace(BigNumber number, Element<Dim>* element1, Element<Dim>* element2, CartesianShape<Dim, Dim - 1>* shape) :
		Poisson_DG_Face<Dim>(number, element1, element2),
		Poisson_HHO_Face<Dim>(number, element1, element2),
		CartesianShape<Dim, Dim - 1>(*shape),
		Face<Dim>(number, element1,  element2)
	{}

public:

	//----------------------------------------------------//
	//                 Face implementation                //
	//----------------------------------------------------//

	void Serialize(ostream& os) const override
	{
		Face<Dim>::Serialize(os);
		os << ": ";
		CartesianShape<Dim, Dim - 1>::Serialize(os);
		/*os << ", finer faces = ";
		for (Face<Dim> * finerFace : this->FinerFaces)
			os << finerFace->Number << ",";*/
	}

	double GetDiameter()
	{
		return max({ CartesianShape<Dim, Dim - 1>::WidthX, CartesianShape<Dim, Dim - 1>::WidthY, CartesianShape<Dim, Dim - 1>::WidthZ });
	}

	double Measure()
	{
		return CartesianShape<Dim, Dim - 1>::Measure;
	}

	RefPoint ConvertToReference(DomPoint domainPoint)
	{
		return CartesianShape<Dim, Dim - 1>::ConvertToReference(domainPoint);
	}
	DomPoint ConvertToDomain(RefPoint referenceElementPoint)
	{
		return CartesianShape<Dim, Dim - 1>::ConvertToDomain(referenceElementPoint);
	}

	double ComputeIntegral(function<double(RefPoint)> func)
	{
		return CartesianShape<Dim, Dim - 1>::ComputeIntegral(func);
	}
	double ComputeIntegral(function<double(RefPoint)> func, int polynomialDegree)
	{
		return CartesianShape<Dim, Dim - 1>::ComputeIntegral(func, polynomialDegree);
	}

	Face<Dim>* CreateSameGeometricFace(BigNumber number, Element<Dim>* element1)
	{
		Face<Dim>* copy = new CartesianFace<Dim>(number, element1, NULL, this);
		copy->IsDomainBoundary = this->IsDomainBoundary;
		return copy;
	}

	void ExportFaceToMatlab(FILE* file)
	{
		fprintf(file, "%llu %.17g %.17g %.17g %.17g %.17g %.17g %d %d\n", this->Number, this->Origin->X, this->Origin->Y, this->Origin->Z, this->WidthX, this->WidthY, this->WidthZ, this->Orientation, this->IsDomainBoundary);
	}

	//----------------------------------------------------------------//
	//                 Poisson_HHO_Face implementation                //
	//----------------------------------------------------------------//

	DenseMatrix FaceMassMatrix(FunctionalBasis<Dim-1>* basis)
	{
		return CartesianShape<Dim, Dim - 1>::FaceMassMatrix(basis);
	}

};