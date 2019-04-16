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

	vector<Point> GetNodalPoints(FunctionalBasis<Dim>* basis)
	{
		vector<Point> points(basis->Size());
		if (Dim == 1)
		{
			GaussLegendre gauss(basis->Size());
			for (int i = 0; i < basis->Size(); ++i)
				points.push_back(gauss.Point(i));
		}
		else if (Dim == 2)
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
		else if (Dim == 3)
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