#pragma once

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
	Point Origin;
	double Width;
	CartesianShapeOrientation Orientation;

	CartesianShape(Point origin, double width) : 
		CartesianShape(origin, width, CartesianShapeOrientation::None)
	{}

	CartesianShape(Point origin, double width, CartesianShapeOrientation orientation)
	{
		this->Origin = origin;
		this->Width = width;
		this->Orientation = orientation;
	}

	double Integral(BasisFunction<ShapeDim>* phi)
	{
		double h = this->Width;
		return pow(h/2, ShapeDim) * Utils::Integral(phi);
	}

	double MassTerm(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		double h = this->Width;

		function<double(Point)> functionToIntegrate = [phi1, phi2](Point p) {
			return phi1->Eval(p)*phi2->Eval(p);
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 2;
		return pow(h / 2, ShapeDim) * Utils::Integral<ShapeDim>(nQuadPoints, functionToIntegrate);
	}

	double IntegralGradGrad(BasisFunction<ShapeDim>* phi1, BasisFunction<ShapeDim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi1->GetDegree() == 0)
			return 0;

		double h = this->Width;

		function<double(Point)> functionToIntegrate = [phi1, phi2](Point p) {
			return Element<ShapeDim>::InnerProduct(phi1->Grad(p), phi2->Grad(p));
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		return pow(h / 2, ShapeDim-2) * Utils::Integral<ShapeDim>(nQuadPoints, functionToIntegrate);
	}

	Point ConvertToDomain(Point referenceElementPoint)
	{
		Point p;
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

	Point ConvertToReference(Point domainPoint)
	{
		Point refPoint;
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
			//double y1 = this->Origin.Y;

			double x2 = x1 + this->Width;
			//double y2 = y1;
			if (this->Orientation == CartesianShapeOrientation::Vertical)
			{
				x2 = x1;
				//y2 = y1 + this->Width;
			}

			double x = domainPoint.X;
			//double y = domainPoint.Y;

			refPoint.X = 2 / (x2 - x1) * x - (x2 + x1) / (x2 - x1);
		}
		else
			assert(false);
		return refPoint;
	}

	vector<Point> GetNodalPoints(FunctionalBasis<ShapeDim>* basis)
	{
		vector<Point> points(basis->Size());
		if (ShapeDim == 1)
		{
			GaussLegendre gauss(basis->Size());
			for (int i = 0; i < basis->Size(); ++i)
				points.push_back(gauss.Point(i));
		}
		else if (ShapeDim == 2)
		{
			if (basis->FullTensorization)
			{
				GaussLegendre gauss(basis->GetDegree() + 1);
				for (int i = 0; i < basis->GetDegree() + 1; ++i)
					for (int j = 0; j < basis->GetDegree() + 1; ++j)
						points.push_back(Point(gauss.Point(i), gauss.Point(j)));
			}
			else
			{
				cout << "Error: Must use the Q space to generate nodal points!" << endl;
				exit(EXIT_FAILURE);
			}
		}
		else if (ShapeDim == 3)
		{
			if (basis->FullTensorization)
			{
				GaussLegendre gauss(basis->GetDegree() + 1);
				for (int i = 0; i < basis->GetDegree() + 1; ++i)
					for (int j = 0; j < basis->GetDegree() + 1; ++j)
						for (int k = 0; k < basis->GetDegree() + 1; ++k)
							points.push_back(Point(gauss.Point(i), gauss.Point(j), gauss.Point(k)));
			}
			else
			{
				cout << "Error: Must use the Q space to generate nodal points!" << endl;
				exit(EXIT_FAILURE);
			}
		}
		return points;
	}

};