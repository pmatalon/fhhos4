#pragma once
#include "ReferenceShape.h"
#include "../TestCases/Diffusion/Tensor.h"
#include <mutex>
using namespace std;

struct GeometricMapping
{
	int NFunctions;
	vector<double> Coeffs;
	vector<double> Exponents;
};

template <int Dim>
class PhysicalShape : public GeometricShape<Dim>
{
private:
	map<RefPoint, DomPoint> _domPoints;
	mutex* _mutex = nullptr;
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

	virtual vector<DomPoint> Vertices() const = 0;

	virtual bool IsConvex() const = 0;

	virtual DomPoint InteriorPoint() const = 0;

	virtual bool Contains(const DomPoint& p) const = 0;

	virtual bool IsDegenerated() const = 0;

	// Transformation to reference element
	virtual DomPoint ConvertToDomain(const RefPoint& refPoint) const = 0; // Mapping
	virtual RefPoint ConvertToReference(const DomPoint& domainPoint) const = 0; // Inverse mapping

	virtual vector<double> MappingCoefficients() const { Utils::FatalError("MappingCoefficients() not implemented for this element"); }

	virtual DimMatrix<Dim> InverseJacobianTranspose(const RefPoint& p) const = 0;
	virtual double DetJacobian(const RefPoint& p) const = 0;
	virtual int DetJacobianDegree() const = 0;

	virtual bool MapsToACartesianShape() const = 0;

	virtual PhysicalShape<Dim>* CreateCopy() const = 0;

	virtual void RefineWithoutCoarseOverlap(const vector<PhysicalShape<Dim - 1>*>& doNotCross) = 0;

	//--------------------//
	//      Geometry      //
	//--------------------//

	virtual bool ConvexHullEmbeds(const PhysicalShape<Dim>* other) const
	{
		if (IsConvex())
		{
			for (const DomPoint& v : other->Vertices())
			{
				if (Contains(v))
					return false;
			}
			return true;
		}
		assert(false);
		return false;
	}

	virtual void ExportToMatlab(string color = "r") const
	{
		assert(false && "To be implemented in the subclass");
	}

	virtual bool IsGeneralPolygon() const
	{
		return false;
	}
	// For squares, triangles, etc, just one subshape: itself.
	// For polyhedra: decomposition into simplices.
	virtual vector<const PhysicalShape<Dim>*> SubShapes() const
	{
		assert(false);
		return {};
	}
	virtual vector<PhysicalShape<Dim>*> SubShapes()
	{
		assert(false);
		return {};
	}

	virtual void ExportSubShapesToMatlab() const
	{
		assert(false);
	}

	virtual vector<const PhysicalShape<Dim>*> RefinedShapes() const
	{
		assert(false);
		return {};
	}
	virtual vector<PhysicalShape<Dim>*> RefinedShapes()
	{
		assert(false);
		return {};
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

private:
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

	virtual double ComputeIntegralKGradGrad(const Tensor<Dim>& K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		RefFunction functionToIntegrate = [this, &K, phi1, phi2](const RefPoint& p) {
			DimMatrix<Dim> invJ = InverseJacobianTranspose(p);
			DimVector<Dim> gradPhi1 = invJ * phi1->Grad(p);
			DimVector<Dim> gradPhi2 = invJ * phi2->Grad(p);
			return (K * gradPhi1).dot(gradPhi2);
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

public:
	virtual double IntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return this->ComputeIntegralGradGrad(phi1, phi2);
	}

	virtual double IntegralKGradGrad(const Tensor<Dim>& K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return this->ComputeIntegralKGradGrad(K, phi1, phi2);
	}

	virtual DenseMatrix MassMatrix(FunctionalBasis<Dim>* basis) const
	{
		return this->ComputeAndReturnMassMatrix(basis);
	}

	DenseMatrix IntegralGradGradMatrix(FunctionalBasis<Dim>* basis) const
	{
		DenseMatrix m(basis->LocalFunctions.size(), basis->LocalFunctions.size());
		for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
			{
				if (phi2->LocalNumber > phi1->LocalNumber)
					break;
				double value = this->IntegralGradGrad(phi1, phi2);
				m(phi1->LocalNumber, phi2->LocalNumber) = value;
				m(phi2->LocalNumber, phi1->LocalNumber) = value;
			}
		}
		return m;
	}

	//----------------------------//
	//             DG             //
	//----------------------------//

	virtual double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const
	{
		return this->ComputeMassTerm(phi1, phi2);
	}

	//-----------------------------//
	//             HHO             //
	//-----------------------------//

	virtual DenseMatrix CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis) const
	{
		return this->ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
	}

	virtual ~PhysicalShape()
	{
		// TODO find a way to deallocate safely
		/*if (_mutex)
		{
			delete _mutex;
			_mutex = nullptr;
		}*/
	}
};