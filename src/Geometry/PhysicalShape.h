#pragma once
#include "ReferenceShape.h"
#include "../Problem/Tensor.h"
#include <mutex>
using namespace std;

template <int Dim>
class PhysicalShape : public GeometricShape<Dim>
{
private:
	map<RefPoint, DomPoint> _domPoints;
	mutex* _mutex;
public:
	PhysicalShape() : GeometricShape<Dim>()
	{
		_mutex = new mutex;
	}

	//-----------------------//
	//   Virtual functions   //
	//-----------------------//

	// Geometry
	virtual ReferenceShape<Dim>* RefShape() const = 0;

	virtual const vector<Vertex*>& Vertices() const = 0;

	virtual bool IsConvex() const = 0;

	virtual bool Contains(const DomPoint& p) const = 0;

	virtual bool IsDegenerated() const = 0;

	// Transformation to reference element
	virtual DomPoint ConvertToDomain(const RefPoint& refPoint) const = 0; // Mapping
	virtual RefPoint ConvertToReference(const DomPoint& domainPoint) const = 0; // Inverse mapping

	virtual DimMatrix<Dim> InverseJacobianTranspose(const RefPoint& p) const = 0;
	virtual double DetJacobian(const RefPoint& p) const = 0;
	virtual int DetJacobianDegree() const = 0;

	virtual PhysicalShape<Dim>* CreateCopy() const = 0;

	//--------------------//
	//      Geometry      //
	//--------------------//

	bool HasVertex(Vertex* v, bool compareCoordinates = false)
	{
		auto vertices = this->Vertices();
		auto it = find_if(vertices.begin(), vertices.end(), [v, compareCoordinates](Vertex* v2) { return v == v2 || (compareCoordinates && *v == *v2); });
		return it != vertices.end();
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

	bool HasOneVertexInCommonWith(PhysicalShape<Dim>* other)
	{
		for (Vertex* v : other->Vertices())
		{
			if (this->HasVertex(v))
				return true;
		}
		return false;
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

	virtual vector<PhysicalShape<Dim>*> IntersectionWith(PhysicalShape<Dim>* other)
	{
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

	//----------------------------------------//
	// Correspondance RefPoint/DomPoint saved //
	//----------------------------------------//

	DomPoint ConvertToDomain(const RefPoint& refPoint, bool getDomPointFromSaved) const
	{
		auto it = _domPoints.find(refPoint);
		if (it != _domPoints.end())
			return it->second;
		return this->ConvertToDomain(refPoint);
	}
	virtual DomPoint ConvertToDomainAndSaveResult(const RefPoint& refPoint, bool lock = false)
	{
		auto it = _domPoints.find(refPoint);
		if (it != _domPoints.end())
			return it->second;
		DomPoint domPoint = this->ConvertToDomain(refPoint);
		if (lock)
			_mutex->lock();
		SaveDomPoint(refPoint, domPoint);
		if (lock)
			_mutex->unlock();
		return domPoint;
	}
	void ComputeAndSaveDomPoint(const RefPoint& refPoint)
	{
		DomPoint domPoint = this->ConvertToDomain(refPoint);
		SaveDomPoint(refPoint, domPoint);
	}
	void ComputeAndSaveQuadraturePoints(int polynomialDegree)
	{
		for (const RefPoint& refPoint : RefShape()->QuadraturePoints(polynomialDegree))
			ConvertToDomainAndSaveResult(refPoint);
	}
	void ComputeAndSaveQuadraturePoints()
	{
		for (const RefPoint& refPoint : RefShape()->QuadraturePoints())
			ConvertToDomainAndSaveResult(refPoint);
	}
private:
	inline void SaveDomPoint(const RefPoint& refPoint, const DomPoint& domPoint)
	{
		_domPoints.insert(pair<RefPoint, DomPoint>(refPoint, domPoint));
	}
public:
	inline void EmptySavedDomPoints()
	{
		_domPoints.clear();
	}

	//-------------------//
	//     Integrals     //
	//-------------------//

	virtual vector<DomPoint> QuadraturePoints() const
	{
		vector<RefPoint> refPoints = RefShape()->QuadraturePoints();
		vector<DomPoint> domPoints;
		for (RefPoint refPoint : refPoints)
			domPoints.push_back(this->ConvertToDomain(refPoint, true));
		return domPoints;
	}

	double Integral(BasisFunction<Dim>* phi) const
	{
		return GeometricShape<Dim>::Integral(phi);
	}

	virtual double Integral(RefFunction f) const override
	{
		RefFunction func = [this, f](const RefPoint& p) {
			return DetJacobian(p) * f(p);
		};
		return RefShape()->Integral(func);
	}

	virtual double Integral(RefFunction f, int polynomialDegree) const override
	{
		RefFunction func = [this, f](const RefPoint& p) {
			return DetJacobian(p) * f(p);
		};
		return RefShape()->Integral(func, polynomialDegree + DetJacobianDegree());
	}

	virtual double Integral(DomFunction globalFunction) const
	{
		RefFunction refFunction = [this, globalFunction](const RefPoint& refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint, true);
			return globalFunction(domainPoint);
		};

		return Integral(refFunction);
	}

	virtual double Integral(DomFunction globalFunction, int polynomialDegree) const
	{
		RefFunction refFunction = [this, globalFunction](const RefPoint& refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint, true);
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

		RefFunction functionToIntegrate = [this, phi1, phi2](const RefPoint& p) {
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

		RefFunction functionToIntegrate = [this, K, phi1, phi2](const RefPoint& p) {
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

	virtual ~PhysicalShape()
	{
		delete _mutex;
	}
};