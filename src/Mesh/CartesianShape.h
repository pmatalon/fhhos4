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
	DomPoint Origin;
	double WidthX = 0;
	double WidthY = 0;
	double WidthZ = 0;
	CartesianShapeOrientation Orientation = CartesianShapeOrientation::None;
	double Measure;
	bool IsRegular; // true if all widths are equal

	static ReferenceCartesianShape<ShapeDim> ReferenceShape;

	//------------------//
	//   Constructors   //
	//------------------//

	CartesianShape(DomPoint origin, double width)
	{
		this->WidthX = ShapeDim >= 1 ? width : 0;
		this->WidthY = ShapeDim >= 2 ? width : 0;
		this->WidthZ = ShapeDim >= 3 ? width : 0;
		this->IsRegular = true;
		Init(origin, CartesianShapeOrientation::None);
	}

	CartesianShape(DomPoint origin, double widthX, double widthY)
	{
		assert(DomainDim == 2 && ShapeDim == 2);
		this->WidthX = widthX;
		this->WidthY = widthY;
		this->IsRegular = (widthX == widthY);
		Init(origin, CartesianShapeOrientation::None);
	}

	CartesianShape(DomPoint origin, double widthX, double widthY, double widthZ)
	{
		assert(DomainDim == 3 && ShapeDim == 3);
		this->WidthX = widthX;
		this->WidthY = widthY;
		this->WidthZ = widthZ;
		this->IsRegular = (widthX == widthY && widthY == widthZ);
		Init(origin, CartesianShapeOrientation::None);
	}

	CartesianShape(DomPoint origin, double width, CartesianShapeOrientation orientation)
	{
		assert(ShapeDim == DomainDim - 1);
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

	CartesianShape(DomPoint origin, double width1, double width2, CartesianShapeOrientation orientation)
	{
		assert(DomainDim == 3 && ShapeDim == DomainDim - 1);
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

	inline void Init(DomPoint origin, CartesianShapeOrientation orientation)
	{
		this->Origin = origin;
		this->Orientation = orientation;
		this->Measure = (this->WidthX != 0 ? this->WidthX : 1) * (this->WidthY != 0 ? this->WidthY : 1) * (this->WidthZ != 0 ? this->WidthZ : 1);
	}

	void Serialize(ostream& os) const
	{
		if (ShapeDim == 2)
			os << (IsRegular ? "square" : "rectangle") << ", ";
		else if (ShapeDim == 3)
			os << (IsRegular ? "cube" : "parallelepipede") << ", ";
		os << "origin = ";
		Origin.Serialize(os, DomainDim);
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

	//---------------//
	//   Integrals   //
	//---------------//

	double IntegralGlobalFunction(function<double(DomPoint)> func)
	{
		if (ShapeDim == 1)
		{
			double x1 = this->Origin.X;
			double x2 = this->Origin.X + this->WidthX;

			return Utils::Integral(func, x1, x2);
		}
		else if (ShapeDim == 2)
		{
			double x1 = this->Origin.X;
			double x2 = this->Origin.X + this->WidthX;
			double y1 = this->Origin.Y;
			double y2 = this->Origin.Y + this->WidthY;

			return Utils::Integral(func, x1, x2, y1, y2);
		}
		else if (ShapeDim == 3)
		{
			double x1 = this->Origin.X;
			double x2 = this->Origin.X + this->WidthX;
			double y1 = this->Origin.Y;
			double y2 = this->Origin.Y + this->WidthY;
			double z1 = this->Origin.Z;
			double z2 = this->Origin.Z + this->WidthZ;

			return Utils::Integral(func, x1, x2, y1, y2, z1, z2);
		}
		else
			assert(false);
	}

	double Integral(BasisFunction<ShapeDim>* phi)
	{
		return RescalingCoeff() * ReferenceShape.ComputeIntegral(phi);
	}

	double ComputeIntegral(function<double(RefPoint)> func)
	{
		double integralOnReferenceShape = ReferenceShape.ComputeIntegral(func);
		return RescalingCoeff() * integralOnReferenceShape;
	}

	double ComputeIntegral(function<double(RefPoint)> func, int polynomialDegree)
	{
		double integralOnReferenceShape = ReferenceShape.ComputeIntegral(func, polynomialDegree);
		return RescalingCoeff() * integralOnReferenceShape;
	}

	double ComputeIntegralGradGrad(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		vector<double> gradTransfo = GradTransformation();

		function<double(RefPoint)> functionToIntegrate = [phi1, phi2, gradTransfo](RefPoint p) {
			vector<double> gradPhi1 = Utils::Multiply<ShapeDim>(gradTransfo, phi1->Grad(p));
			vector<double> gradPhi2 = Utils::Multiply<ShapeDim>(gradTransfo, phi2->Grad(p));
			return Element<ShapeDim>::InnerProduct(gradPhi1, gradPhi2);
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return ComputeIntegral(functionToIntegrate, polynomialDegree);
	}

	//--------//
	//   DG   //
	//--------//

	double MassTerm(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		return RescalingCoeff() * ReferenceShape.MassTerm(phi1, phi2);
	}

	double StiffnessTerm(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		if (this->IsRegular)
		{
			vector<double> gradTransfo = GradTransformation();
			return RescalingCoeff() * pow(gradTransfo[0], 2) * ReferenceShape.StiffnessTerm(phi1, phi2);
		}
		else
			return ComputeIntegralGradGrad(phi1, phi2);
	}

	//---------//
	//   HHO   //
	//---------//

	Eigen::MatrixXd FaceMassMatrix(FunctionalBasis<ShapeDim>* basis)
	{
		return RescalingCoeff() * ReferenceShape.FaceMassMatrix(basis);
	}

	Eigen::MatrixXd CellMassMatrix(FunctionalBasis<ShapeDim>* basis)
	{
		return RescalingCoeff() * ReferenceShape.CellMassMatrix(basis);
	}

	Eigen::MatrixXd CellReconstructMassMatrix(FunctionalBasis<ShapeDim>* cellBasis, FunctionalBasis<ShapeDim>* reconstructBasis)
	{
		return RescalingCoeff() * ReferenceShape.CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	double IntegralGradGradReconstruct(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		if (this->IsRegular)
		{
			vector<double> gradTransfo = GradTransformation();
			return RescalingCoeff() * pow(gradTransfo[0], 2) * ReferenceShape.ReconstructStiffnessTerm(phi1, phi2);
		}
		else
			return ComputeIntegralGradGrad(phi1, phi2);
	}

private:
	inline double RescalingCoeff()
	{
		return this->Measure / pow(2, ShapeDim);
	}

public:
	DomPoint ConvertToDomain(RefPoint referenceElementPoint)
	{
		DomPoint p;
		if (ShapeDim == DomainDim)
		{
			double x1 = this->Origin.X;
			double x2 = this->Origin.X + this->WidthX;
			double y1 = this->Origin.Y;
			double y2 = this->Origin.Y + this->WidthY;
			double z1 = this->Origin.Z;
			double z2 = this->Origin.Z + this->WidthZ;

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
			p = this->Origin;
		}
		else if (ShapeDim == 1 && DomainDim == 2)
		{
			double x1 = this->Origin.X;
			double y1 = this->Origin.Y;

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
			double x1 = this->Origin.X;
			double y1 = this->Origin.Y;
			double z1 = this->Origin.Z;

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

	RefPoint ConvertToReference(DomPoint domainPoint)
	{
		RefPoint refPoint;
		if (ShapeDim == DomainDim)
		{
			double x1 = this->Origin.X;
			double x2 = this->Origin.X + this->WidthX;
			double y1 = this->Origin.Y;
			double y2 = this->Origin.Y + this->WidthY;
			double z1 = this->Origin.Z;
			double z2 = this->Origin.Z + this->WidthZ;

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
			/*double x1 = this->Origin.X;
			double x2 = x1 + this->WidthX;
			if (this->Orientation == CartesianShapeOrientation::Vertical)
				x2 = x1;

			double x = domainPoint.X;

			refPoint.X = 2 / (x2 - x1) * x - (x2 + x1) / (x2 - x1);*/
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

	vector<double> GradTransformation()
	{
		vector<double> gradTransfo(ShapeDim);
		if (ShapeDim == DomainDim)
		{
			if (ShapeDim >= 1)
				gradTransfo[0] = 2 / this->WidthX;
			if (ShapeDim >= 2)
				gradTransfo[1] = 2 / this->WidthY;
			if (ShapeDim == 3)
				gradTransfo[2] = 2 / this->WidthZ;
		}
		else if (ShapeDim == 1 && DomainDim == 2)
		{
			if (this->Orientation == CartesianShapeOrientation::Horizontal)
				gradTransfo[0] = 2 / this->WidthX;
			else if (this->Orientation == CartesianShapeOrientation::Vertical)
				gradTransfo[0] = 2 / this->WidthY;
			else
				assert(false);
		}
		else if (ShapeDim == 2 && DomainDim == 3)
		{
			assert(false && "GradTransformation: 3D case to be implemented!");
		}
		else
			assert(false);
		return gradTransfo;
	}

	vector<RefPoint> GetNodalPoints(FunctionalBasis<ShapeDim>* basis)
	{
		vector<RefPoint> points;
		points.reserve(basis->Size());
		if (ShapeDim == 1)
		{
			if (basis->GetDegree() == 0)
				points.push_back(RefPoint(0));
			else
			{
				double h = 2 / ((double)basis->Size() - 1);
				for (int i = 0; i < basis->Size(); ++i)
					points.push_back(RefPoint(-1 + i * h));
			}
		}
		else if (ShapeDim == 2)
		{
			if (basis->GetDegree() == 0)
				points.push_back(RefPoint(0, 0));
			else if (basis->FullTensorization)
			{
				double h = 2 / (double)basis->GetDegree();
				for (int i = 0; i < basis->GetDegree() + 1; ++i)
					for (int j = 0; j < basis->GetDegree() + 1; ++j)
						points.push_back(RefPoint(-1 + i * h, -1 + j * h));
			}
			else
			{
				if (basis->GetDegree() == 1)
				{
					points.push_back(RefPoint(-1, 1)); // top-left corner
					points.push_back(RefPoint(1, 1)); // top-right corner
					points.push_back(RefPoint(0, -1)); // bottom-middle
				}
				else if (basis->GetDegree() == 2)
				{
					points.push_back(RefPoint(-1, 1)); // top-left corner
					points.push_back(RefPoint(1, 1)); // top-right corner
					points.push_back(RefPoint(1, -1)); // bottom-right corner
					points.push_back(RefPoint(-1, -1)); // bottom-left corner

					points.push_back(RefPoint(-1, 0)); // center-left
					points.push_back(RefPoint(1, 0)); // center-right
				}
				else if (basis->GetDegree() == 3)
				{
					points.push_back(RefPoint(-1, 1)); // top-left corner
					points.push_back(RefPoint(1, 1)); // top-right corner
					points.push_back(RefPoint(1, -1)); // bottom-right corner
					points.push_back(RefPoint(-1, -1)); // bottom-left corner

					points.push_back(RefPoint(-1, 0)); // center-left
					points.push_back(RefPoint(1, 0)); // center-right
					points.push_back(RefPoint(0, 1)); // top-middle
					points.push_back(RefPoint(0, -1)); // bottom-middle

					points.push_back(RefPoint(-0.5, 0));
					points.push_back(RefPoint(0.5, 0));
				}
				else
				{
					cout << "Warning: interpolation points to be implemented!" << endl;
					int nPoint1D = ceil(sqrt(basis->Size()));
					points.reserve(pow(nPoint1D, 2));
					GaussLegendre gauss(nPoint1D);
					for (int i = 0; i < nPoint1D; ++i)
						for (int j = 0; j < nPoint1D; ++j)
							points.push_back(RefPoint(gauss.Point(i), gauss.Point(j)));
				}
			}
		}
		else if (ShapeDim == 3)
		{
			cout << "Warning: interpolation points to be implemented!" << endl;
			int nPoint1D = ceil(cbrt(basis->Size()));
			points.reserve(pow(nPoint1D, 3));
			GaussLegendre gauss(nPoint1D);
			for (int i = 0; i < nPoint1D; ++i)
				for (int j = 0; j < nPoint1D; ++j)
					for (int k = 0; k < nPoint1D; ++k)
						points.push_back(RefPoint(gauss.Point(i), gauss.Point(j), gauss.Point(k)));
		}
		return points;
	}

};

template <int DomainDim, int ShapeDim>
ReferenceCartesianShape<ShapeDim> CartesianShape<DomainDim, ShapeDim>::ReferenceShape = ReferenceCartesianShape<ShapeDim>();