#pragma once
#include "../BiHarmonicTestCase.h"
using namespace std;

class CubeBiHarTestCase : public BiHarmonicTestCase<3>
{
private:
	ProblemArguments _pb;
public:
	CubeBiHarTestCase(ProblemArguments pb) :
		BiHarmonicTestCase(),
		_pb(pb)
	{
		if (pb.SourceCode.compare("") == 0)
		{
			pb.SourceCode = "poly";
			_pb.SourceCode = "poly";
		}

		// Boundary conditions
		if (pb.SourceCode.compare("exp") == 0)
		{
			this->DirichletBC.Type = PbBoundaryConditions::FullDirichlet;
			this->DirichletBC.BoundaryConditionPartition = BoundaryConditions::DirichletEverywhere;
			this->DirichletBC.DirichletFunction = this->ExpSolution_Dirichlet;
			this->DirichletBC.NeumannFunction = nullptr;

			this->NeumannBC.Type = PbBoundaryConditions::FullNeumann;
			this->NeumannBC.BoundaryConditionPartition = BoundaryConditions::NeumannEverywhere;
			this->NeumannBC.NeumannFunction = this->ExpSolution_Neumann;
			this->NeumannBC.DirichletFunction = nullptr;
		}
		else
		{
			this->DirichletBC = BoundaryConditions::HomogeneousDirichletEverywhere();
			this->NeumannBC = BoundaryConditions::HomogeneousNeumannEverywhere();
		}

		// Source function
		if (pb.SourceCode.compare("poly") == 0)
			this->SourceFunction = this->PolySource;
		else if (pb.SourceCode.compare("exp") == 0)
			this->SourceFunction = this->ExpSource;
		else if (pb.SourceCode.compare("one") == 0)
			this->SourceFunction = Utils::ConstantFunctionOne;
		else
			Utils::FatalError("Unmanaged source code (use 'poly', 'exp' or 'one')");

		// Exact solution
		if (pb.SourceCode.compare("poly") == 0)
		{
			this->ExactSolution = this->PolySolution;
			this->MinusLaplacianOfSolution = this->PolyMinusLaplacianOfSolution;
			this->MinusLaplacianOfSolution_Dirichlet = this->PolyMinusLaplacianOfSolution_Dirichlet;
		}
		else if (pb.SourceCode.compare("exp") == 0)
		{
			this->ExactSolution = this->ExpSolution;
			this->MinusLaplacianOfSolution = this->ExpMinusLaplacianOfSolution;
			this->MinusLaplacianOfSolution_Dirichlet = this->ExpMinusLaplacianOfSolution_Dirichlet;
		}
	}
	/*
	//--------------//
	//     Sine     //
	//--------------//
	static double SineSolution(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return pow(sin(M_PI * x), 2) * pow(sin(M_PI * y), 2);
	}
	static double SineMinusLaplacianOfSolution(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 4 * pow(M_PI, 2) * pow(sin(M_PI * x), 2) * pow(sin(M_PI * y), 2) - 2 * pow(M_PI, 2) * pow(cos(M_PI * x), 2) * pow(sin(M_PI * y), 2) - 2 * pow(M_PI, 2) * pow(cos(M_PI * y), 2) * pow(sin(M_PI * x), 2);
	}
	static double SineMinusLaplacianOfSolution_Dirichlet(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		if (abs(x) < Utils::NumericalZero || abs(x - 1) < Utils::NumericalZero) // x = 0 || x = 1
			return -2 * pow(M_PI, 2) * pow(sin(M_PI * y), 2);
		else if (abs(y) < Utils::NumericalZero || abs(y - 1) < Utils::NumericalZero) // y = 0 || y = 1
			return -2 * pow(M_PI, 2) * pow(sin(M_PI * x), 2);
		Utils::FatalError("Should never happen");
		return 0;
	}
	static double SineSource(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		return 24 * pow(M_PI, 4) * pow(sin(M_PI * x), 2) * pow(sin(M_PI * y),2) + 8 * pow(M_PI, 4) * pow(cos(M_PI * x),2) * pow(cos(M_PI * y),2) - 16 * pow(M_PI, 4) * pow(cos(M_PI * x),2) * pow(sin(M_PI * y),2) - 16 * pow(M_PI, 4) * pow(cos(M_PI * y),2) * pow(sin(M_PI * x),2);
	}*/

	//--------------//
	//     Poly     //
	//--------------//
	static double PolySolution(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return pow(x, 4) * pow(x - 1, 2) * pow(y, 4) * pow(y - 1, 2) * pow(z, 4) * pow(z - 1, 2);
	}
	static double PolyMinusLaplacianOfSolution(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return -2 * pow(x,4) * pow(y,4) * pow(z,4) * pow(x - 1, 2) * pow(y - 1, 2) - 2 * pow(x,4) * pow(y,4) * pow(z,4) * pow(x - 1, 2) * pow(z - 1, 2) - 2 * pow(x,4) * pow(y,4) * pow(z,4) * pow(y - 1, 2) * pow(z - 1, 2) - 8 * pow(x,3) * pow(y,4) * pow(z,4) * (2 * x - 2) * pow(y - 1, 2) * pow(z - 1, 2) - 8 * pow(x,4) * pow(y,3) * pow(z,4) * (2 * y - 2) * pow(x - 1, 2) * pow(z - 1, 2) - 8 * pow(x,4) * pow(y,4) * pow(z,3) * (2 * z - 2) * pow(x - 1, 2) * pow(y - 1, 2) - 12 * pow(x,2) * pow(y,4) * pow(z,4) * pow(x - 1, 2) * pow(y - 1, 2) * pow(z - 1, 2) - 12 * pow(x,4) * pow(y,2) * pow(z,4) * pow(x - 1, 2) * pow(y - 1, 2) * pow(z - 1, 2) - 12 * pow(x,4) * pow(y,4) * pow(z,2) * pow(x - 1, 2) * pow(y - 1, 2) * pow(z - 1, 2);
	}
	static double PolyMinusLaplacianOfSolution_Dirichlet(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		if (     abs(x)     < Utils::NumericalZero) // x = 0 (left)
			return 0;
		else if (abs(x - 1) < Utils::NumericalZero) // x = 1 (right)
			return (y * y * y * y) * (z * z * z * z) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * -2.0;
		else if (abs(y)     < Utils::NumericalZero) // y = 0 (front)
			return 0;
		else if (abs(y - 1) < Utils::NumericalZero) // y = 1 (back)
			return (x * x * x * x) * (z * z * z * z) * pow(x - 1.0, 2.0) * pow(z - 1.0, 2.0) * -2.0;
		else if (abs(z)     < Utils::NumericalZero) // z = 0 (bottom)
			return 0;
		else if (abs(z - 1) < Utils::NumericalZero) // z = 1 (top)
			return (x * x * x * x) * (y * y * y * y) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * -2.0;
		Utils::FatalError("Should never happen");
		return 0;
	}
	static double PolySource(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return (x * x * x * x) * (y * y * y * y) * (z * z * z * z) * pow(x - 1.0, 2.0) * 8.0 + (x * x * x * x) * (y * y * y * y) * (z * z * z * z) * pow(y - 1.0, 2.0) * 8.0 + (x * x * x * x) * (y * y * y * y) * (z * z * z * z) * pow(z - 1.0, 2.0) * 8.0 + (x * x * x) * (y * y * y * y) * (z * z * z * z) * (x * 2.0 - 2.0) * pow(y - 1.0, 2.0) * 3.2E1 + (x * x * x * x) * (y * y * y) * (z * z * z * z) * (y * 2.0 - 2.0) * pow(x - 1.0, 2.0) * 3.2E1 + (x * x) * (y * y * y * y) * (z * z * z * z) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * 4.8E1 + (x * x * x * x) * (y * y) * (z * z * z * z) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * 4.8E1 + (x * x * x * x) * (y * y * y * y) * (z * z) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * 1.44E2 + (x * x * x) * (y * y * y * y) * (z * z * z * z) * (x * 2.0 - 2.0) * pow(z - 1.0, 2.0) * 3.2E1 + (x * x * x * x) * (y * y * y * y) * (z * z * z) * (z * 2.0 - 2.0) * pow(x - 1.0, 2.0) * 3.2E1 + (x * x) * (y * y * y * y) * (z * z * z * z) * pow(x - 1.0, 2.0) * pow(z - 1.0, 2.0) * 4.8E1 + (x * x * x * x) * (y * y) * (z * z * z * z) * pow(x - 1.0, 2.0) * pow(z - 1.0, 2.0) * 1.44E2 + (x * x * x * x) * (y * y * y * y) * (z * z) * pow(x - 1.0, 2.0) * pow(z - 1.0, 2.0) * 4.8E1 + (x * x * x * x) * (y * y * y) * (z * z * z * z) * (y * 2.0 - 2.0) * pow(z - 1.0, 2.0) * 3.2E1 + (x * x * x * x) * (y * y * y * y) * (z * z * z) * (z * 2.0 - 2.0) * pow(y - 1.0, 2.0) * 3.2E1 + (x * x) * (y * y * y * y) * (z * z * z * z) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * 1.44E2 + (x * x * x * x) * (y * y) * (z * z * z * z) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * 4.8E1 + (x * x * x * x) * (y * y * y * y) * (z * z) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * 4.8E1 + (x * x * x * x) * (y * y * y * y) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * 2.4E1 + (x * x * x * x) * (z * z * z * z) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * 2.4E1 + (y * y * y * y) * (z * z * z * z) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * 2.4E1 + x * (y * y * y * y) * (z * z * z * z) * (x * 2.0 - 2.0) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * 9.6E1 + (x * x * x * x) * y * (z * z * z * z) * (y * 2.0 - 2.0) * pow(x - 1.0, 2.0) * pow(z - 1.0, 2.0) * 9.6E1 + (x * x * x * x) * (y * y * y * y) * z * (z * 2.0 - 2.0) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * 9.6E1 + (x * x * x) * (y * y * y) * (z * z * z * z) * (x * 2.0 - 2.0) * (y * 2.0 - 2.0) * pow(z - 1.0, 2.0) * 1.28E2 + (x * x * x) * (y * y * y * y) * (z * z * z) * (x * 2.0 - 2.0) * (z * 2.0 - 2.0) * pow(y - 1.0, 2.0) * 1.28E2 + (x * x * x * x) * (y * y * y) * (z * z * z) * (y * 2.0 - 2.0) * (z * 2.0 - 2.0) * pow(x - 1.0, 2.0) * 1.28E2 + (x * x) * (y * y * y) * (z * z * z * z) * (y * 2.0 - 2.0) * pow(x - 1.0, 2.0) * pow(z - 1.0, 2.0) * 1.92E2 + (x * x) * (y * y * y * y) * (z * z * z) * (z * 2.0 - 2.0) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * 1.92E2 + (x * x * x) * (y * y) * (z * z * z * z) * (x * 2.0 - 2.0) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * 1.92E2 + (x * x * x) * (y * y * y * y) * (z * z) * (x * 2.0 - 2.0) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * 1.92E2 + (x * x * x * x) * (y * y) * (z * z * z) * (z * 2.0 - 2.0) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * 1.92E2 + (x * x * x * x) * (y * y * y) * (z * z) * (y * 2.0 - 2.0) * pow(x - 1.0, 2.0) * pow(z - 1.0, 2.0) * 1.92E2 + (x * x) * (y * y) * (z * z * z * z) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * 2.88E2 + (x * x) * (y * y * y * y) * (z * z) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * 2.88E2 + (x * x * x * x) * (y * y) * (z * z) * pow(x - 1.0, 2.0) * pow(y - 1.0, 2.0) * pow(z - 1.0, 2.0) * 2.88E2;
	}

	//-------------//
	//     Exp     //
	//-------------//
	static double ExpSolution(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return x * z * sin(M_PI * y) * exp(-x * y);
	}
	static double ExpSolution_Dirichlet(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		if (     abs(x)     < Utils::NumericalZero) // x = 0 (left)
			return 0;
		else if (abs(x - 1) < Utils::NumericalZero) // x = 1 (right)
			return z * exp(-y) * sin(M_PI * y);
		else if (abs(y)     < Utils::NumericalZero) // y = 0 (front)
			return 0;
		else if (abs(y - 1) < Utils::NumericalZero) // y = 1 (back)
			return 0;
		else if (abs(z)     < Utils::NumericalZero) // z = 0 (bottom)
			return 0;
		else if (abs(z - 1) < Utils::NumericalZero) // z = 1 (top)
			return x * exp(-x * y) * sin(M_PI * y);
		Utils::FatalError("Should never happen");
		return 0;
	}
	static double ExpSolution_Neumann(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		if (     abs(x)     < Utils::NumericalZero) // x = 0 (left)
			return -z * sin(M_PI * y);
		else if (abs(x - 1) < Utils::NumericalZero) // x = 1 (right)
			return exp(-y) * z * sin(M_PI * y) - exp(-y) * y * z * sin(M_PI * y);
		else if (abs(y)     < Utils::NumericalZero) // y = 0 (front)
			return -M_PI * x * z;
		else if (abs(y - 1) < Utils::NumericalZero) // y = 1 (back)
			return -M_PI * exp(-x) * x * z;
		else if (abs(z)     < Utils::NumericalZero) // z = 0 (bottom)
			return -exp(-x * y) * x * sin(M_PI * y);
		else if (abs(z - 1) < Utils::NumericalZero) // z = 1 (top)
			return exp(-x * y) * x * sin(M_PI * y);
		Utils::FatalError("Should never happen");
		return 0;
	}
	static double ExpMinusLaplacianOfSolution(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return 2 * y * z * exp(-x * y) * sin(M_PI * y) - pow(x,3) * z * exp(-x * y) * sin(M_PI * y) + 2 * pow(x,2) * z * M_PI * exp(-x * y) * cos(M_PI * y) + x * z * pow(M_PI,2) * exp(-x * y) * sin(M_PI * y) - x * pow(y,2) * z * exp(-x * y) * sin(M_PI * y);
	}
	static double ExpMinusLaplacianOfSolution_Dirichlet(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		if (         abs(x) < Utils::NumericalZero) // x = 0 (left)
			return 2 * y * z * sin(M_PI * y);
		else if (abs(x - 1) < Utils::NumericalZero) // x = 1 (right)
			return 2 * y * z * exp(-y) * sin(M_PI * y) - z * exp(-y) * sin(M_PI * y) + z * pow(M_PI,2) * exp(-y) * sin(M_PI * y) - pow(y,2) * z * exp(-y) * sin(M_PI * y) + 2 * z * M_PI * exp(-y) * cos(M_PI * y);
		else if (abs(y)     < Utils::NumericalZero) // y = 0 (front)
			return 2 * M_PI * pow(x,2) * z;
		else if (abs(y - 1) < Utils::NumericalZero) // y = 1 (back)
			return -2 * pow(x,2) * z * M_PI * exp(-x);
		else if (abs(z)     < Utils::NumericalZero) // z = 0 (bottom)
			return 0;
		else if (abs(z - 1) < Utils::NumericalZero) // z = 1 (top)
			return 2 * y * exp(-x * y) * sin(M_PI * y) - pow(x,3) * exp(-x * y) * sin(M_PI * y) - x * pow(y,2) * exp(-x * y) * sin(M_PI * y) + 2 * pow(x,2) * M_PI * exp(-x * y) * cos(M_PI * y) + x * pow(M_PI,2) * exp(-x * y) * sin(M_PI * y);
		Utils::FatalError("Should never happen");
		return 0;
	}
	static double ExpSource(const DomPoint& p)
	{
		double x = p.X;
		double y = p.Y;
		double z = p.Z;
		return pow(x,5) * z * exp(-x * y) * sin(M_PI * y) - 4 * pow(y,3) * z * exp(-x * y) * sin(M_PI * y) - 8 * z * M_PI * exp(-x * y) * cos(M_PI * y) + 12 * x * z * exp(-x * y) * sin(M_PI * y) - 4 * pow(x,4) * z * M_PI * exp(-x * y) * cos(M_PI * y) + x * z * pow(M_PI,4) * exp(-x * y) * sin(M_PI * y) + 4 * y * z * pow(M_PI,2) * exp(-x * y) * sin(M_PI * y) - 12 * pow(x,2) * y * z * exp(-x * y) * sin(M_PI * y) + x * pow(y,4) * z * exp(-x * y) * sin(M_PI * y) + 4 * pow(x,2) * z * pow(M_PI,3) * exp(-x * y) * cos(M_PI * y) - 6 * pow(x,3) * z * pow(M_PI,2) * exp(-x * y) * sin(M_PI * y) + 2 * pow(x,3) * pow(y,2) * z * exp(-x * y) * sin(M_PI * y) - 4 * pow(x,2) * pow(y,2) * z * M_PI * exp(-x * y) * cos(M_PI * y) - 2 * x * pow(y,2) * z * pow(M_PI,2) * exp(-x * y) * sin(M_PI * y) + 16 * x * y * z * M_PI * exp(-x * y) * cos(M_PI * y);
	}





	string Code() override
	{
		if (_pb.SourceCode.compare("sine") == 0)
			return "sine";
		else if (_pb.SourceCode.compare("poly") == 0)
			return "poly";
		else if (_pb.SourceCode.compare("exp") == 0)
			return "exp";
		else if (_pb.SourceCode.compare("one") == 0)
			return "one";
		else
			return _pb.SourceCode;
	}
	string Description() override
	{
		if (_pb.SourceCode.compare("sine") == 0)
			return "sine solution";
		else if (_pb.SourceCode.compare("poly") == 0)
			return "polynomial solution";
		else if (_pb.SourceCode.compare("exp") == 0)
			return "exponential solution";
		else if (_pb.SourceCode.compare("one") == 0)
			return "one";
		else
			return _pb.SourceCode;
	}
};
