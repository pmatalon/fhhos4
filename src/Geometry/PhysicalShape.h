#pragma once
#include "ReferenceShape.h"
#include "../Problem/Tensor.h"

template <int Dim>
class PhysicalShape : public GeometricShape<Dim>
{
public:
	PhysicalShape() : GeometricShape<Dim>() {}

	//-----------------------//
	//   Virtual functions   //
	//-----------------------//

	// Geometry
	virtual ReferenceShape<Dim>* RefShape() const = 0;

	virtual vector<Vertex*> Vertices() const = 0;

	virtual bool IsConvex() const = 0;

	virtual bool Contains(DomPoint p) const = 0;

	virtual bool IsDegenerated() const = 0;

	// Transformation to reference element
	virtual DomPoint ConvertToDomain(RefPoint refPoint) const = 0; // Mapping
	virtual RefPoint ConvertToReference(DomPoint domainPoint) const = 0; // Inverse mapping

	virtual DimMatrix<Dim> InverseJacobianTranspose(RefPoint p) const = 0;
	virtual double DetJacobian(RefPoint p) const = 0;
	virtual int DetJacobianDegree() const = 0;

	virtual PhysicalShape<Dim>* CreateCopy() const = 0;

	//--------------------//
	//      Geometry      //
	//--------------------//

	bool HasVertex(Vertex* v, bool compareCoordinates = false)
	{
		for (Vertex* v2 : this->Vertices())
		{
			if (v == v2 || (compareCoordinates && *v == *v2))
				return true;
		}
		return false;
	}

	bool HasSameVertices(PhysicalShape<Dim>* other, bool compareCoordinates = false)
	{
		if (this->Vertices().size() != other->Vertices().size())
			return false;

		for (Vertex* v : this->Vertices())
		{
			if (!other->HasVertex(v, compareCoordinates))
				return false;
		}
		return true;
	}

	bool Contains(PhysicalShape<Dim>* other)
	{
		if (!this->Contains(other->Center()))
			return false;

		for (auto v : other->Vertices())
		{
			if (!this->HasVertex(v, true) && !this->Contains(*v))
				return false;
		}

		return true;
	}

	virtual bool ConvexHullEmbeds(PhysicalShape<Dim>* other) const
	{
		if (IsConvex())
		{
			for (Vertex* v : other->Vertices())
			{
				if (Contains(*v))
					return false;
			}
			return true;
		}
		assert(false);
	}

	virtual void ExportToMatlab(string color = "r") const
	{
		assert(false && "To be implemented in the subclass");
	}

	// For squares, triangles, etc, just one subshape: itself.
	// For polyhedra: decomposition into simplices.
	// For agglomerates: list of shapes forming the agglomerate.
	virtual bool IsMadeOfSubShapes() const
	{
		return false;
	}
	virtual vector<PhysicalShape<Dim>*> SubShapes() const
	{
		assert(false);
	}
	virtual void ExportSubShapesToMatlab() const
	{
		assert(false);
	}

	virtual void ReshapeByMovingIntersection(Vertex* oldIntersect, Vertex* newIntersect)
	{
		assert(false && "To be implemented in subclasses");
	}

	bool IsIn(const vector<PhysicalShape<Dim>*>& list)
	{
		for (PhysicalShape<Dim>* s : list)
		{
			if (s == this)
				return true;
		}
		return false;
	}

	//-------------------//
	//     Integrals     //
	//-------------------//

	virtual vector<DomPoint> QuadraturePoints() const
	{
		vector<RefPoint> refPoints = RefShape()->QuadraturePoints();
		vector<DomPoint> domPoints;
		for (RefPoint refPoint : refPoints)
			domPoints.push_back(this->ConvertToDomain(refPoint));
		return domPoints;
	}

	double Integral(BasisFunction<Dim>* phi) const
	{
		return GeometricShape<Dim>::Integral(phi);
	}

	virtual double Integral(RefFunction f) const override
	{
		RefFunction func = [this, f](RefPoint p) {
			return DetJacobian(p) * f(p);
		};
		return RefShape()->Integral(func);
	}

	virtual double Integral(RefFunction f, int polynomialDegree) const override
	{
		RefFunction func = [this, f](RefPoint p) {
			return DetJacobian(p) * f(p);
		};
		return RefShape()->Integral(func, polynomialDegree + DetJacobianDegree());
	}

	virtual double Integral(DomFunction globalFunction) const
	{
		RefFunction refFunction = [this, globalFunction](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return globalFunction(domainPoint);
		};

		return Integral(refFunction);
	}

	virtual double Integral(DomFunction globalFunction, int polynomialDegree) const
	{
		RefFunction refFunction = [this, globalFunction](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return globalFunction(domainPoint);
		};

		return Integral(refFunction, polynomialDegree);
	}

	//----------------------------//
	//     Specific integrals     //
	//----------------------------//

	virtual double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		RefFunction functionToIntegrate = [this, phi1, phi2](RefPoint p) {
			DimMatrix<Dim> invJ = InverseJacobianTranspose(p);
			DimVector<Dim> gradPhi1 = invJ * phi1->Grad(p);
			DimVector<Dim> gradPhi2 = invJ * phi2->Grad(p);
			return gradPhi1.dot(gradPhi2);
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

	virtual double ComputeIntegralKGradGrad(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		RefFunction functionToIntegrate = [this, K, phi1, phi2](RefPoint p) {
			DimMatrix<Dim> invJ = InverseJacobianTranspose(p);
			DimVector<Dim> gradPhi1 = invJ * phi1->Grad(p);
			DimVector<Dim> gradPhi2 = invJ * phi2->Grad(p);
			return (K * gradPhi1).dot(gradPhi2);
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

	//----------------------------//
	//             DG             //
	//----------------------------//

	virtual double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->ComputeMassTerm(phi1, phi2);
	}

	virtual double StiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->ComputeIntegralGradGrad(phi1, phi2);
	}

	//-----------------------------//
	//             HHO             //
	//-----------------------------//

	virtual DenseMatrix FaceMassMatrix(FunctionalBasis<Dim>* basis)
	{
		return this->ComputeAndReturnMassMatrix(basis);
	}

	virtual DenseMatrix CellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		return this->ComputeAndReturnMassMatrix(basis);
	}

	virtual DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		return this->ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
	}

	virtual double IntegralKGradGradReconstruct(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->ComputeIntegralKGradGrad(K, phi1, phi2);
	}
};