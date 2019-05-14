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
	double Width;
	CartesianShapeOrientation Orientation;

	static ReferenceCartesianShape<ShapeDim> ReferenceShape;

	CartesianShape(DomPoint origin, double width) : 
		CartesianShape(origin, width, CartesianShapeOrientation::None)
	{}

	CartesianShape(DomPoint origin, double width, CartesianShapeOrientation orientation)
	{
		this->Origin = origin;
		this->Width = width;
		this->Orientation = orientation;
	}

	double Measure()
	{
		return pow(this->Width, ShapeDim);
	}

	double IntegralGlobalFunction(function<double(DomPoint)> func)
	{
		if (ShapeDim == 1)
		{
			double x1 = this->Origin.X;
			double x2 = this->Origin.X + this->Width;

			return Utils::Integral(func, x1, x2);
		}
		else if (ShapeDim == 2)
		{
			double x1 = this->Origin.X;
			double x2 = this->Origin.X + this->Width;
			double y1 = this->Origin.Y;
			double y2 = this->Origin.Y + this->Width;

			return Utils::Integral(func, x1, x2, y1, y2);
		}
		else if (ShapeDim == 3)
		{
			double x1 = this->Origin.X;
			double x2 = this->Origin.X + this->Width;
			double y1 = this->Origin.Y;
			double y2 = this->Origin.Y + this->Width;
			double z1 = this->Origin.Z;
			double z2 = this->Origin.Z + this->Width;

			return Utils::Integral(func, x1, x2, y1, y2, z1, z2);
		}
		else
			assert(false);
	}

	double Integral(BasisFunction<ShapeDim>* phi)
	{
		return Rescale(ReferenceShape.ComputeIntegral(phi), 0);
	}

	double ComputeIntegral(function<double(RefPoint)> func, int numberOfDerivatives)
	{
		double integralOnReferenceShape = ReferenceShape.ComputeIntegral(func);
		return Rescale(integralOnReferenceShape, numberOfDerivatives);
	}

	double ComputeIntegral(function<double(RefPoint)> func, int numberOfDerivatives, int polynomialDegree)
	{
		double integralOnReferenceShape = ReferenceShape.ComputeIntegral(func, polynomialDegree);
		return Rescale(integralOnReferenceShape, numberOfDerivatives);
	}

	Eigen::MatrixXd ComputeAndReturnFaceMassMatrix(FunctionalBasis<ShapeDim>* basis)
	{
		return RescaleMass(ReferenceShape.StoredFaceMassMatrix());
	}

	Eigen::MatrixXd ComputeAndReturnCellMassMatrix(FunctionalBasis<ShapeDim>* basis)
	{
		return RescaleMass(ReferenceShape.StoredCellMassMatrix());
	}

	Eigen::MatrixXd ComputeAndReturnCellReconstructMassMatrix(FunctionalBasis<ShapeDim>* cellBasis, FunctionalBasis<ShapeDim>* reconstructBasis)
	{
		return RescaleMass(ReferenceShape.StoredCellReconstructMassMatrix());
	}


	double MassTerm(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		return RescaleMass(ReferenceShape.MassTerm(phi1, phi2));
	}

	double IntegralGradGrad(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi1->GetDegree() == 0)
			return 0;

		return RescaleStiffness(ReferenceShape.StiffnessTerm(phi1, phi2));
	}

	double IntegralGradGradReconstruct(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi1->GetDegree() == 0)
			return 0;

		return RescaleStiffness(ReferenceShape.ReconstructStiffnessTerm(phi1, phi2));
	}

	double ComputeIntegralGradGrad(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi1->GetDegree() == 0)
			return 0;

		return RescaleStiffness(ReferenceShape.ComputeIntegralGradGrad(phi1, phi2));
	}

private:
	Eigen::MatrixXd Rescale(const Eigen::MatrixXd& matrixOnReferenceElement, int numberOfDerivatives)
	{
		double h = this->Width;
		return pow(h / 2, ShapeDim - numberOfDerivatives) * matrixOnReferenceElement;
	}

	double Rescale(double termOnReferenceElement, int numberOfDerivatives)
	{
		double h = this->Width;
		return pow(h / 2, ShapeDim - numberOfDerivatives) * termOnReferenceElement;
	}

	Eigen::MatrixXd RescaleMass(const Eigen::MatrixXd& massMatrixOnReferenceElement)
	{
		return Rescale(massMatrixOnReferenceElement, 0);
	}

	double RescaleMass(double massTermOnReferenceElement)
	{
		return Rescale(massTermOnReferenceElement, 0);
	}

	Eigen::MatrixXd RescaleStiffness(const Eigen::MatrixXd& stiffnessMatrixOnReferenceElement)
	{
		return Rescale(stiffnessMatrixOnReferenceElement, 2);
	}

	double RescaleStiffness(double stiffnessTermOnReferenceElement)
	{
		return Rescale(stiffnessTermOnReferenceElement, 2);
	}

public:
	DomPoint ConvertToDomain(RefPoint referenceElementPoint)
	{
		DomPoint p;
		if (ShapeDim == DomainDim)
		{
			double x1 = this->Origin.X;
			double x2 = this->Origin.X + this->Width;
			double y1 = this->Origin.Y;
			double y2 = this->Origin.Y + this->Width;
			double z1 = this->Origin.Z;
			double z2 = this->Origin.Z + this->Width;

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
		else if (ShapeDim == 1 && DomainDim == 2)
		{
			double x1 = this->Origin.X;
			double y1 = this->Origin.Y;

			double x2 = x1 + this->Width;
			double y2 = y1;
			if (this->Orientation == CartesianShapeOrientation::Vertical)
			{
				x2 = x1;
				y2 = y1 + this->Width;
			}

			double t = referenceElementPoint.X;
			p.X = (x2 - x1) / 2 * t + (x2 + x1) / 2;
			p.Y = (y2 - y1) / 2 * t + (y2 + y1) / 2;
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
			double x2 = this->Origin.X + this->Width;
			double y1 = this->Origin.Y;
			double y2 = this->Origin.Y + this->Width;
			double z1 = this->Origin.Z;
			double z2 = this->Origin.Z + this->Width;

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
			double x1 = this->Origin.X;
			double x2 = x1 + this->Width;
			if (this->Orientation == CartesianShapeOrientation::Vertical)
				x2 = x1;

			double x = domainPoint.X;

			refPoint.X = 2 / (x2 - x1) * x - (x2 + x1) / (x2 - x1);
		}
		else
			assert(false);
		return refPoint;
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
				cout << "Warning: these interpolation points have never been tested." << endl;
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