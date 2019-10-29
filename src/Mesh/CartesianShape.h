#pragma once
#include "ReferenceCartesianShape.h"

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
class CartesianShape
{
public:
	DomPoint* Origin;
	double WidthX = 0;
	double WidthY = 0;
	double WidthZ = 0;
	CartesianShapeOrientation Orientation = CartesianShapeOrientation::None;
	double Measure;
	DomPoint Center;
	bool IsRegular; // true if all widths are equal

	static ReferenceCartesianShape<ShapeDim> ReferenceShape;

	//------------------//
	//   Constructors   //
	//------------------//

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
		this->Measure = (this->WidthX != 0 ? this->WidthX : 1) * (this->WidthY != 0 ? this->WidthY : 1) * (this->WidthZ != 0 ? this->WidthZ : 1);
		if (ShapeDim == 1)
		{
			if (DomainDim == 1)
				this->Center = DomPoint(Origin->X + WidthX / 2);
			else if (this->Orientation == CartesianShapeOrientation::Horizontal)
				this->Center = DomPoint(Origin->X + WidthX / 2, Origin->Y);
			else
				this->Center = DomPoint(Origin->X, Origin->Y + WidthY / 2);
		}
		else if (ShapeDim == 2)
		{
			if (DomainDim == 2)
			{
				if (this->Orientation == CartesianShapeOrientation::Horizontal)
					this->Center = DomPoint(Origin->X + WidthX / 2, Origin->Y);
				else
					this->Center = DomPoint(Origin->X, Origin->Y + WidthY / 2);
			}
			else
				assert(false && "not implemented yet");
		}
		else if (ShapeDim == 3)
		{
			if (this->Orientation == CartesianShapeOrientation::InXOY)
				this->Center = DomPoint(Origin->X + WidthX / 2, Origin->Y + WidthY / 2, Origin->Z);
			else if (this->Orientation == CartesianShapeOrientation::InXOZ)
				this->Center = DomPoint(Origin->X + WidthX / 2, Origin->Y, Origin->Z + WidthZ / 2);
			else if (this->Orientation == CartesianShapeOrientation::InYOZ)
				this->Center = DomPoint(Origin->X, Origin->Y + WidthY / 2, Origin->Z + WidthZ / 2);
		}
	}

	void Serialize(ostream& os) const
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

	inline double DetJacobian() const
	{
		return this->Measure / ReferenceShape.Measure();
	}

	DimMatrix<ShapeDim> InverseJacobianTranspose() const
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

	DomPoint ConvertToDomain(RefPoint referenceElementPoint) const
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

	RefPoint ConvertToReference(DomPoint domainPoint) const
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

	//-------------------------------------//
	//              Integrals              //
	//-------------------------------------//

	double Integral(DomFunction func) const
	{
		if (ShapeDim == 1)
		{
			double x1 = this->Origin->X;
			double x2 = this->Origin->X + this->WidthX;

			return Utils::Integral(func, x1, x2);
		}
		else if (ShapeDim == 2)
		{
			double x1 = this->Origin->X;
			double x2 = this->Origin->X + this->WidthX;
			double y1 = this->Origin->Y;
			double y2 = this->Origin->Y + this->WidthY;

			return Utils::Integral(func, x1, x2, y1, y2);
		}
		else if (ShapeDim == 3)
		{
			double x1 = this->Origin->X;
			double x2 = this->Origin->X + this->WidthX;
			double y1 = this->Origin->Y;
			double y2 = this->Origin->Y + this->WidthY;
			double z1 = this->Origin->Z;
			double z2 = this->Origin->Z + this->WidthZ;

			return Utils::Integral(func, x1, x2, y1, y2, z1, z2);
		}
		else
			assert(false);
	}

	double Integral(BasisFunction<ShapeDim>* phi) const
	{
		return DetJacobian() * ReferenceShape.Integral(phi);
	}
	double Integral(RefFunction func) const
	{
		return DetJacobian() * ReferenceShape.Integral(func);
	}
	double Integral(RefFunction func, int polynomialDegree) const
	{
		return DetJacobian() * ReferenceShape.Integral(func, polynomialDegree);
	}

	double ComputeIntegralGradGrad(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2) const
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		DimMatrix<ShapeDim> invJ = InverseJacobianTranspose();

		RefFunction functionToIntegrate = [phi1, phi2, invJ](RefPoint p) {
			DimVector<ShapeDim> gradPhi1 = invJ * phi1->Grad(p);
			DimVector<ShapeDim> gradPhi2 = invJ * phi2->Grad(p);
			return gradPhi1.dot(gradPhi2);
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

	double ComputeIntegralKGradGrad(Tensor<ShapeDim>* K, BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2) const
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		DimMatrix<ShapeDim> invJ = InverseJacobianTranspose();

		RefFunction functionToIntegrate = [K, phi1, phi2, invJ](RefPoint p) {
			DimVector<ShapeDim> gradPhi1 = invJ * phi1->Grad(p);
			DimVector<ShapeDim> gradPhi2 = invJ * phi2->Grad(p);
			return (K * gradPhi1).dot(gradPhi2);
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

	//----------------------------//
	//             DG             //
	//----------------------------//

	double MassTerm(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		return DetJacobian() * ReferenceShape.MassTerm(phi1, phi2);
	}

	double StiffnessTerm(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		if (this->IsRegular)
		{
			DimMatrix<ShapeDim> invJ = InverseJacobianTranspose();
			return DetJacobian() * pow(invJ(0, 0), 2) * ReferenceShape.StiffnessTerm(phi1, phi2);
		}
		else
			return ComputeIntegralGradGrad(phi1, phi2);
	}

	//-----------------------------//
	//             HHO             //
	//-----------------------------//

	DenseMatrix FaceMassMatrix(FunctionalBasis<ShapeDim>* basis)
	{
		return DetJacobian() * ReferenceShape.FaceMassMatrix(basis);
	}

	DenseMatrix CellMassMatrix(FunctionalBasis<ShapeDim>* basis)
	{
		return DetJacobian() * ReferenceShape.CellMassMatrix(basis);
	}

	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<ShapeDim>* cellBasis, FunctionalBasis<ShapeDim>* reconstructBasis)
	{
		return DetJacobian() * ReferenceShape.CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	double IntegralKGradGradReconstruct(Tensor<ShapeDim>* K, BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		if (this->IsRegular)
		{
			DimMatrix<ShapeDim> invJ = InverseJacobianTranspose();
			return DetJacobian() * pow(invJ(0, 0), 2) * ReferenceShape.ReconstructKStiffnessTerm(K, phi1, phi2);
		}
		else
			return ComputeIntegralKGradGrad(K, phi1, phi2);
	}
};

template <int DomainDim, int ShapeDim>
ReferenceCartesianShape<ShapeDim> CartesianShape<DomainDim, ShapeDim>::ReferenceShape = ReferenceCartesianShape<ShapeDim>();