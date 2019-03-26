#pragma once

template <int Dim>
class CartesianShape
{
public:
	Point Origin;
	double Width;

	CartesianShape(Point origin, double width)
	{
		this->Origin = origin;
		this->Width = width;
	}

	double Integral(BasisFunction<Dim>* phi)
	{
		double h = this->Width;
		return pow(h/2, Dim) * Utils::Integral(phi);
	}

	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		double h = this->Width;

		function<double(Point)> functionToIntegrate = [phi1, phi2](Point p) {
			return phi1->Eval(p)*phi2->Eval(p);
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree() + 2;
		return pow(h / 2, Dim) * Utils::Integral<Dim>(nQuadPoints, functionToIntegrate);
	}

	double IntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi1->GetDegree() == 0)
			return 0;

		double h = this->Width;

		function<double(Point)> functionToIntegrate = [phi1, phi2](Point p) {
			return Element<Dim>::InnerProduct(phi1->Grad(p), phi2->Grad(p));
		};

		int nQuadPoints = phi1->GetDegree() + phi2->GetDegree();
		return pow(h / 2, Dim-2) * Utils::Integral<Dim>(nQuadPoints, functionToIntegrate);
	}

	Point ConvertToDomain(Point referenceElementPoint)
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

		Point p;
		p.X = (x2 - x1) / 2 * t + (x2 + x1) / 2;
		if (Dim >= 2)
		{
			p.Y = (y2 - y1) / 2 * u + (y2 + y1) / 2;
		}
		if (Dim == 3)
		{
			p.Z = (z2 - z1) / 2 * v + (z2 + z1) / 2;
		}
		return p;
	}

};