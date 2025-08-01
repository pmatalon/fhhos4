#pragma once
#include "ReferenceCartesianShape.h"
#include "PhysicalShapeWithConstantJacobian.h"

enum class CartesianShapeOrientation : unsigned
{
	None = 0,
	Horizontal = 1,
	Vertical = 2,
	InXOY = 3,
	InXOZ = 4,
	InYOZ = 5
};

template <int DomainDim, int ShapeDim = DomainDim>
class CartesianShape : public PhysicalShapeWithConstantJacobian<ShapeDim>
{
private:
	vector<DomPoint> _vertices;

	double _diameter;
	double _measure;
	DomPoint _center;
	double _inRadius = 0;
public:
	DomPoint* Origin;
	double WidthX = 0;
	double WidthY = 0;
	double WidthZ = 0;
	CartesianShapeOrientation Orientation = CartesianShapeOrientation::None;
	bool IsRegular; // true if all widths are equal

	static ReferenceCartesianShape<ShapeDim> RefCartShape;

	//------------------//
	//   Constructors   //
	//------------------//

	CartesianShape() {}

	CartesianShape(const CartesianShape &) = default; // Copy constructor

	CartesianShape(DomPoint* origin, double width)
	{
		assert(width > 0);
		this->WidthX = ShapeDim >= 1 ? width : 0;
		this->WidthY = ShapeDim >= 2 ? width : 0;
		this->WidthZ = ShapeDim >= 3 ? width : 0;
		this->IsRegular = true;
		Init(origin, CartesianShapeOrientation::None);
	}

	CartesianShape(DomPoint* origin, double widthX, double widthY)
	{
		assert(DomainDim == 2 && ShapeDim == 2);
		assert(widthX > 0 && widthY > 0);
		this->WidthX = widthX;
		this->WidthY = widthY;
		this->IsRegular = (widthX == widthY);
		Init(origin, CartesianShapeOrientation::None);
	}

	CartesianShape(DomPoint* origin, double widthX, double widthY, double widthZ)
	{
		assert(DomainDim == 3 && ShapeDim == 3);
		assert(widthX > 0 && widthY > 0 && widthZ > 0);
		this->WidthX = widthX;
		this->WidthY = widthY;
		this->WidthZ = widthZ;
		this->IsRegular = (widthX == widthY && widthY == widthZ);
		Init(origin, CartesianShapeOrientation::None);
	}

	CartesianShape(DomPoint* origin, double width, CartesianShapeOrientation orientation)
	{
		assert(ShapeDim == DomainDim - 1);
		assert(ShapeDim == 0 || width > 0);
		if (DomainDim == 2)
		{
			if (orientation == CartesianShapeOrientation::Horizontal)
				this->WidthX = width;
			else if (orientation == CartesianShapeOrientation::Vertical)
				this->WidthY = width;
			else
				assert(false);
		}
		else if (DomainDim == 3)
		{
			if (orientation == CartesianShapeOrientation::InXOY)
			{
				this->WidthX = width;
				this->WidthY = width;
			}
			else if (orientation == CartesianShapeOrientation::InXOZ)
			{
				this->WidthX = width;
				this->WidthZ = width;
			}
			else if (orientation == CartesianShapeOrientation::InYOZ)
			{
				this->WidthY = width;
				this->WidthZ = width;
			}
			else
				assert(false);
		}
		this->IsRegular = true;
		Init(origin, orientation);
	}

	CartesianShape(DomPoint* origin, double width1, double width2, CartesianShapeOrientation orientation)
	{
		assert(DomainDim == 3 && ShapeDim == DomainDim - 1);
		assert(width1 > 0 && width2 > 0);
		if (orientation == CartesianShapeOrientation::InXOY)
		{
			this->WidthX = width1;
			this->WidthY = width2;
		}
		else if (orientation == CartesianShapeOrientation::InXOZ)
		{
			this->WidthX = width1;
			this->WidthZ = width2;
		}
		else if (orientation == CartesianShapeOrientation::InYOZ)
		{
			this->WidthY = width1;
			this->WidthZ = width2;
		}
		else
			assert(false);
		this->IsRegular = false;
		Init(origin, orientation);
	}

	inline void Init(DomPoint* origin, CartesianShapeOrientation orientation)
	{
		this->Origin = origin;
		this->Orientation = orientation;
		this->_diameter = max({ WidthX, WidthY, WidthZ });
		this->_measure = (this->WidthX != 0 ? this->WidthX : 1) * (this->WidthY != 0 ? this->WidthY : 1) * (this->WidthZ != 0 ? this->WidthZ : 1);
		if (ShapeDim == 1)
		{
			if (DomainDim == 1)
				this->_center = DomPoint(Origin->X + WidthX / 2);
			else if (this->Orientation == CartesianShapeOrientation::Horizontal)
				this->_center = DomPoint(Origin->X + WidthX / 2, Origin->Y);
			else
				this->_center = DomPoint(Origin->X, Origin->Y + WidthY / 2);
		}
		else if (ShapeDim == 2)
		{
			if (DomainDim == 2)
			{
				this->_center = DomPoint(Origin->X + WidthX / 2, Origin->Y + WidthY / 2);
				_inRadius = min(WidthX, WidthY);
			}
			else
			{
				if (this->Orientation == CartesianShapeOrientation::InXOY)
					this->_center = DomPoint(Origin->X + WidthX / 2, Origin->Y + WidthY / 2, Origin->Z);
				else if (this->Orientation == CartesianShapeOrientation::InXOZ)
					this->_center = DomPoint(Origin->X + WidthX / 2, Origin->Y, Origin->Z + WidthZ / 2);
				else if (this->Orientation == CartesianShapeOrientation::InYOZ)
					this->_center = DomPoint(Origin->X, Origin->Y + WidthY / 2, Origin->Z + WidthZ / 2);
			}
		}
		else if (ShapeDim == 3)
		{
			assert(DomainDim == 3);
			this->_center = DomPoint(Origin->X + WidthX / 2, Origin->Y + WidthY / 2, Origin->Z + WidthZ / 2);
			_inRadius = min({ WidthX, WidthY, WidthZ });
		}
	}

	//-------------------------------------//
	//       Geometric information         //
	//-------------------------------------//

	PhysicalShape<ShapeDim>* CreateCopy() const override
	{
		return new CartesianShape<DomainDim,ShapeDim>(*this);
	}

	ReferenceShape<ShapeDim>* RefShape() const override
	{
		return &RefCartShape;
	}

	bool MapsToACartesianShape() const override
	{
		return true;
	}

	// Not a great solution, the vertices should pass by the constructor,
	// but hey, a little "quick and dirty" can't do no harm bro...!
	inline void SetVertices(const vector<DomPoint>& vertices)
	{
		_vertices = vertices;
	}

	inline vector<DomPoint> Vertices() const override
	{
		return _vertices;
	}

	bool IsDegenerated() const override
	{
		Utils::FatalError("The function CartesianShape::IsDegenerated() is not implemented.");
		return true;
	}

	static ReferenceCartesianShape<ShapeDim>* InitReferenceShape()
	{
		return &RefCartShape;
	}

	double Diameter() const override
	{
		return _diameter;
	}
	double Measure() const override
	{
		return _measure;
	}
	DomPoint Center() const override
	{
		return _center;
	}
	inline DomPoint InteriorPoint() const override
	{
		return _center;
	}
	inline bool IsConvex() const override
	{
		return true;
	}
	double InRadius() const override
	{
		return _inRadius;
	}
	inline bool Contains(const DomPoint& p) const override
	{
		if (DomainDim == 1 && ShapeDim == 1)
			return (Origin->X <= p.X && p.X <= Origin->X + WidthX);
		else if (DomainDim == 2 && ShapeDim == 2)
			return (Origin->X <= p.X && p.X <= Origin->X + WidthX) && (Origin->Y <= p.Y && p.Y <= Origin->Y + WidthY);
		else if (DomainDim == 3 && ShapeDim == 3)
			return (Origin->X <= p.X && p.X <= Origin->X + WidthX) && (Origin->Y <= p.Y && p.Y <= Origin->Y + WidthY) && (Origin->Z <= p.Z && p.Z <= Origin->Z + WidthZ);
		else
		{
			Utils::FatalError("Contains() not implemented");
			return false;
		}
	}

	void RefineWithoutCoarseOverlap(const vector<PhysicalShape<ShapeDim - 1>*>& doNotCross) override
	{
		Utils::FatalError("TO BE IMPLEMENTED");
	}

	void Serialize(ostream& os) const override
	{
		if (ShapeDim == 2)
			os << (IsRegular ? "square" : "rectangle") << ", ";
		else if (ShapeDim == 3)
			os << (IsRegular ? "cube" : "parallelepipede") << ", ";
		os << "origin = ";
		Origin->Serialize(os, DomainDim);
		if (ShapeDim == DomainDim)
		{
			if (IsRegular)
				os << ", width = " << WidthX;
			else
			{
				os << ", widthX = " << WidthX;
				if (DomainDim >= 2)
					os << ", widthY = " << WidthY;
				else if (DomainDim >= 3)
					os << ", widthZ = " << WidthZ;
			}
		}
		else
		{
			if (DomainDim == 1)
				os << ", width = " << WidthX;
			else if (DomainDim == 2)
			{
				if (WidthX != 0)
					os << ", widthX = " << WidthX;
				else
					os << ", widthY = " << WidthY;
			}
			else if (DomainDim == 3)
			{
				if (WidthX != 0)
					os << ", widthX = " << WidthX;
				if (WidthY != 0)
					os << ", widthY = " << WidthY;
				if (WidthZ != 0)
					os << ", widthZ = " << WidthZ;
			}
		}
	}

	//-------------------------------------//
	// Transformation to reference element //
	//-------------------------------------//

	inline double DetJacobian() const override
	{
		return _measure / RefCartShape.Measure();
	}

	DimMatrix<ShapeDim> InverseJacobianTranspose() const override
	{
		DimMatrix<ShapeDim> invJ = DimMatrix<ShapeDim>::Zero();
		if (ShapeDim == DomainDim)
		{
			if (ShapeDim >= 1)
				invJ(0, 0) = 2 / this->WidthX;
			if (ShapeDim >= 2)
				invJ(1, 1) = 2 / this->WidthY;
			if (ShapeDim == 3)
				invJ(2, 2) = 2 / this->WidthZ;
		}
		else if (ShapeDim == 1 && DomainDim == 2)
		{
			if (this->Orientation == CartesianShapeOrientation::Horizontal)
				invJ(0, 0) = 2 / this->WidthX;
			else if (this->Orientation == CartesianShapeOrientation::Vertical)
				invJ(0, 0) = 2 / this->WidthY;
			else
				assert(false);
		}
		else if (ShapeDim == 2 && DomainDim == 3)
		{
			assert(false && "InverseJacobianTranspose: 3D case to be implemented!");
		}
		else
			assert(false);
		return invJ;
	}

	DomPoint ConvertToDomain(const RefPoint& referenceElementPoint) const override
	{
		DomPoint p;
		if (ShapeDim == DomainDim)
		{
			double x1 = this->Origin->X;
			double x2 = this->Origin->X + this->WidthX;
			double y1 = this->Origin->Y;
			double y2 = this->Origin->Y + this->WidthY;
			double z1 = this->Origin->Z;
			double z2 = this->Origin->Z + this->WidthZ;

			double t = referenceElementPoint.X;
			double u = referenceElementPoint.Y;
			double v = referenceElementPoint.Z;

			if (ShapeDim >= 1)
				p.X = (x2 - x1) / 2 * t + (x2 + x1) / 2;
			if (ShapeDim >= 2)
				p.Y = (y2 - y1) / 2 * u + (y2 + y1) / 2;
			if (ShapeDim == 3)
				p.Z = (z2 - z1) / 2 * v + (z2 + z1) / 2;
		}
		else if (ShapeDim == 0 && DomainDim == 1)
		{
			p = *(this->Origin);
		}
		else if (ShapeDim == 1 && DomainDim == 2)
		{
			double x1 = this->Origin->X;
			double y1 = this->Origin->Y;

			double x2 = x1 + this->WidthX;
			double y2 = y1;
			if (this->Orientation == CartesianShapeOrientation::Vertical)
			{
				x2 = x1;
				y2 = y1 + this->WidthY;
			}

			double t = referenceElementPoint.X;
			p.X = (x2 - x1) / 2 * t + (x2 + x1) / 2;
			p.Y = (y2 - y1) / 2 * t + (y2 + y1) / 2;
		}
		else if (ShapeDim == 2 && DomainDim == 3)
		{
			double x1 = this->Origin->X;
			double y1 = this->Origin->Y;
			double z1 = this->Origin->Z;

			double x2 = x1 + this->WidthX;
			double y2 = y1 + this->WidthY;
			double z2 = z1 + this->WidthZ;

			double t = referenceElementPoint.X;
			double u = referenceElementPoint.Y;

			if (this->Orientation == CartesianShapeOrientation::InXOY)
			{
				p.X = (x2 - x1) / 2 * t + (x2 + x1) / 2;
				p.Y = (y2 - y1) / 2 * u + (y2 + y1) / 2;
				p.Z = z1;
			}
			else if (this->Orientation == CartesianShapeOrientation::InXOZ)
			{
				p.X = (x2 - x1) / 2 * t + (x2 + x1) / 2;
				p.Y = y1;
				p.Z = (z2 - z1) / 2 * u + (z2 + z1) / 2;
			}
			else if (this->Orientation == CartesianShapeOrientation::InYOZ)
			{
				p.X = x1;
				p.Y = (y2 - y1) / 2 * t + (y2 + y1) / 2;
				p.Z = (z2 - z1) / 2 * u + (z2 + z1) / 2;
			}
			else
				assert(false);
		}
		else
			assert(false);
		return p;
	}

	RefPoint ConvertToReference(const DomPoint& domainPoint) const override
	{
		RefPoint refPoint;
		if (ShapeDim == DomainDim)
		{
			double x1 = this->Origin->X;
			double x2 = this->Origin->X + this->WidthX;
			double y1 = this->Origin->Y;
			double y2 = this->Origin->Y + this->WidthY;
			double z1 = this->Origin->Z;
			double z2 = this->Origin->Z + this->WidthZ;

			double x = domainPoint.X;
			double y = domainPoint.Y;
			double z = domainPoint.Z;

			if (ShapeDim >= 1)
				refPoint.X = 2 / (x2 - x1) * x - (x2 + x1) / (x2 - x1);
			if (ShapeDim >= 2)
				refPoint.Y = 2 / (y2 - y1) * y - (y2 + y1) / (y2 - y1);
			if (ShapeDim == 3)
				refPoint.Z = 2 / (z2 - z1) * z - (z2 + z1) / (z2 - z1);
		}
		else if (ShapeDim == 1 && DomainDim == 2)
		{
			if (this->Orientation == CartesianShapeOrientation::Horizontal)
			{
				double x1 = this->Origin->X;
				double x2 = x1 + this->WidthX;
				double x = domainPoint.X;
				refPoint.X = 2 / (x2 - x1) * x - (x2 + x1) / (x2 - x1);
			}
			else if (this->Orientation == CartesianShapeOrientation::Vertical)
			{
				double y1 = this->Origin->Y;
				double y2 = y1 + this->WidthY;
				double y = domainPoint.Y;
				refPoint.X = 2 / (y2 - y1) * y - (y2 + y1) / (y2 - y1);
			}
			else
				assert(false);
		}
		else if (ShapeDim == 2 && DomainDim == 3)
		{
			assert(false && "ConvertToReference: 3D case to be implemented!");
		}
		else
			assert(false);
		return refPoint;
	}

	static GeometricMapping MappingInfo()
	{
		GeometricMapping mapping;
		mapping.NFunctions = 3; // 1, t, u

		mapping.Coeffs = vector<double>(mapping.NFunctions * mapping.NFunctions);
		for (int i = 0; i < mapping.NFunctions; i++)
			mapping.Coeffs[i * mapping.NFunctions + i] = 1; // Identity

		//                    t,u,v
		mapping.Exponents = { 0,0,0,   // 1   = t^0 * u^0 * v^0
							  1,0,0,   // t   = t^1 * u^0 * v^0
							  0,1,0 }; // u   = t^0 * u^1 * v^0
		return mapping;
	}

	vector<double> MappingCoefficients() const override
	{
		// refer to ConvertToDomain(). The basis functions are 1, t, u.
		double x1 = this->Origin->X;
		double y1 = this->Origin->Y;
		double z1 = this->Origin->Z;

		double x2 = x1 + this->WidthX;
		double y2 = y1 + this->WidthY;
		double z2 = z1 + this->WidthZ;

		if (ShapeDim == 2 && DomainDim == 3)
		{
			if (this->Orientation == CartesianShapeOrientation::InXOY)
			{
				//                   1,             t,             u
				return { (x2 + x1) / 2, (x2 - x1) / 2,             0,   // X
						 (y2 + y1) / 2,             0, (y2 - y1) / 2,   // Y
						            z1,             0,             0 }; // Z
			}
			else if (this->Orientation == CartesianShapeOrientation::InXOZ)
			{
				//                   1,             t,             u
				return { (x2 + x1) / 2, (x2 - x1) / 2,             0,   // X
						            y1,             0,             0,   // Y
						 (z2 + z1) / 2,             0, (z2 - z1) / 2 }; // Z
			}
			else if (this->Orientation == CartesianShapeOrientation::InYOZ)
			{
				//                   1,             t,             u
				return {            x1,             0,             0,   // X
					  	 (y2 + y1) / 2, (y2 - y1) / 2,             0,   // Y
						 (z2 + z1) / 2,             0, (z2 - z1) / 2 }; // Z
			}
		}
		Utils::FatalError("MappingCoefficients() not implemented for CartesianShape<" + to_string(DomainDim) + "," + to_string(ShapeDim) + ">.");
		return vector<double>();
	}

	//-------------------------------------//
	//              Integrals              //
	//-------------------------------------//

	double Integral(DomFunction func) const override
	{
		if (ShapeDim == 1)
		{
			double x1 = this->Origin->X;
			double x2 = this->Origin->X + this->WidthX;

			return Integral(func, x1, x2);
		}
		else if (ShapeDim == 2)
		{
			double x1 = this->Origin->X;
			double x2 = this->Origin->X + this->WidthX;
			double y1 = this->Origin->Y;
			double y2 = this->Origin->Y + this->WidthY;

			return Integral(func, x1, x2, y1, y2);
		}
		else if (ShapeDim == 3)
		{
			double x1 = this->Origin->X;
			double x2 = this->Origin->X + this->WidthX;
			double y1 = this->Origin->Y;
			double y2 = this->Origin->Y + this->WidthY;
			double z1 = this->Origin->Z;
			double z2 = this->Origin->Z + this->WidthZ;

			return Integral(func, x1, x2, y1, y2, z1, z2);
		}
		else
		{
			assert(false);
			return 0.0;
		}
	}

	//----------------------------//
	//             DG             //
	//----------------------------//

	double IntegralGradGrad(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2) const override
	{
		if (this->IsRegular)
		{
			DimMatrix<ShapeDim> invJ = InverseJacobianTranspose();
			return DetJacobian() * pow(invJ(0, 0), 2) * RefCartShape.StiffnessTerm(phi1, phi2);
		}
		else
			return PhysicalShape<ShapeDim>::IntegralGradGrad(phi1, phi2);
	}

	//-----------------------------//
	//             HHO             //
	//-----------------------------//
	
	double IntegralKGradGrad(const Tensor<ShapeDim>& K, BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2) const override
	{
		if (this->IsRegular)
		{
			DimMatrix<ShapeDim> invJ = InverseJacobianTranspose();
			return DetJacobian() * pow(invJ(0, 0), 2) * RefCartShape.ReconstructStiffnessTerm(K, phi1, phi2);
		}
		else
			return PhysicalShape<ShapeDim>::IntegralKGradGrad(K, phi1, phi2);
	}

private:

	//-------------//
	// Integral 1D //
	//-------------//

	static double Integral(int nPoints, std::function<double(double)> func, double x1, double x2)
	{
		GaussLegendre* gs = GaussLegendre::Get(nPoints);
		return gs->Quadrature(func, x1, x2);
	}

	inline static double Integral(std::function<double(double)> func, double x1, double x2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2);
	}

	static double Integral(int nPoints, DomFunction func, double x1, double x2)
	{
		function<double(double)> funcToIntegrate = [func](double x) {
			return func(DomPoint(x));
		};
		return Integral(nPoints, funcToIntegrate, x1, x2);
	}

	inline static double Integral(DomFunction func, double x1, double x2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2);
	}

	//-------------//
	// Integral 2D //
	//-------------//

	// Integral on [x1, x2] x [y1, y2]
	static double Integral(int nPoints, std::function<double(double, double)> func, double x1, double x2, double y1, double y2)
	{
		GaussLegendre* gs = GaussLegendre::Get(nPoints);
		return gs->Quadrature(func, x1, x2, y1, y2);
	}

	inline static double Integral(std::function<double(double, double)> func, double x1, double x2, double y1, double y2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2, y1, y2);
	}

	static double Integral(int nPoints, DomFunction func, double x1, double x2, double y1, double y2)
	{
		function<double(double, double)> funcToIntegrate = [func](double x, double y) {
			return func(DomPoint(x, y));
		};
		return Integral(nPoints, funcToIntegrate, x1, x2, y1, y2);
	}

	inline static double Integral(DomFunction func, double x1, double x2, double y1, double y2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2, y1, y2);
	}

	//-------------//
	// Integral 3D //
	//-------------//

	// Integral on [x1, x2] x [y1, y2] x [z1, z2]
	static double Integral(int nPoints, std::function<double(double, double, double)> func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		GaussLegendre* gs = GaussLegendre::Get(nPoints);
		return gs->Quadrature(func, x1, x2, y1, y2, z1, z2);
	}

	inline static double Integral(std::function<double(double, double, double)> func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2, y1, y2, z1, z2);
	}

	static double Integral(int nPoints, DomFunction func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		function<double(double, double, double)> funcToIntegrate = [func](double x, double y, double z) {
			return func(DomPoint(x, y, z));
		};
		return Integral(nPoints, funcToIntegrate, x1, x2, y1, y2, z1, z2);
	}

	inline static double Integral(DomFunction func, double x1, double x2, double y1, double y2, double z1, double z2)
	{
		return Integral(GaussLegendre::MAX_POINTS, func, x1, x2, y1, y2, z1, z2);
	}

	//---------------------------------------------------------------------//
	// This is f***ing useless, it should be automatic due to inheritance! //
	// But without that it doesn't compile for some reason :-(             //
	//---------------------------------------------------------------------//
public:
	double Integral(RefFunction func) const override
	{
		return PhysicalShapeWithConstantJacobian<ShapeDim>::Integral(func);
	}
	double Integral(RefFunction func, int polynomialDegree) const override
	{
		return PhysicalShapeWithConstantJacobian<ShapeDim>::Integral(func, polynomialDegree);
	}
};

template <int DomainDim, int ShapeDim>
ReferenceCartesianShape<ShapeDim> CartesianShape<DomainDim, ShapeDim>::RefCartShape = ReferenceCartesianShape<ShapeDim>();

using RectangleShape = CartesianShape<2>;