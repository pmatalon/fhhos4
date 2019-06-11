#pragma once
#include "Element.h"
#include "CartesianShape.h"
#include "../HHO/Poisson_HHO_Element.h"

template <int Dim>
class CartesianElement : public Poisson_DG_Element<Dim>, public Poisson_HHO_Element<Dim>, public CartesianShape<Dim>
{
public:

	CartesianElement(BigNumber number, DomPoint origin, double width) :
		Poisson_DG_Element<Dim>(number),
		Poisson_HHO_Element<Dim>(number), 
		CartesianShape<Dim>(origin, width)
	{}

	CartesianElement(BigNumber number, DomPoint origin, double widthX, double widthY) :
		Poisson_DG_Element<Dim>(number),
		Poisson_HHO_Element<Dim>(number),
		CartesianShape<Dim>(origin, widthX, widthY)
	{}

	CartesianElement(BigNumber number, DomPoint origin, double widthX, double widthY, double widthZ) :
		Poisson_DG_Element<Dim>(number),
		Poisson_HHO_Element<Dim>(number),
		CartesianShape<Dim>(origin, widthX, widthY, widthZ)
	{}

	void Serialize(ostream& os) const override
	{
		Element<Dim>::Serialize(os);
		os << ", ";
		CartesianShape<Dim>::Serialize(os);
	}

	double GetDiameter() override
	{
		return max({ CartesianShape<Dim>::WidthX, CartesianShape<Dim>::WidthY, CartesianShape<Dim>::WidthZ });
	}

	double Measure()
	{
		return CartesianShape<Dim>::Measure;
	}

	double DiffusionCoefficient(DiffusionPartition diffusionPartition) override
	{
		DomPoint origin = CartesianShape<Dim>::Origin;
		return diffusionPartition.Coefficient(origin);
	}

	double IntegralGlobalFunction(function<double(DomPoint)> func)
	{
		return CartesianShape<Dim>::IntegralGlobalFunction(func);
	}

	double Integral(BasisFunction<Dim>* phi)
	{
		return CartesianShape<Dim>::Integral(phi);
	}

	double ComputeIntegral(function<double(RefPoint)> func)
	{
		return CartesianShape<Dim>::ComputeIntegral(func);
	}

	double ComputeIntegral(function<double(RefPoint)> func, int polynomialDegree)
	{
		return CartesianShape<Dim>::ComputeIntegral(func, polynomialDegree);
	}

	double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return CartesianShape<Dim>::ComputeIntegralGradGrad(phi1, phi2);
	}

	DomPoint ConvertToDomain(RefPoint referenceElementPoint)
	{
		return CartesianShape<Dim>::ConvertToDomain(referenceElementPoint);
	}

	RefPoint ConvertToReference(DomPoint domainPoint)
	{
		return CartesianShape<Dim>::ConvertToReference(domainPoint);
	}

	vector<double> GradTransformation()
	{
		return CartesianShape<Dim>::GradTransformation();
	}

	vector<RefPoint> GetNodalPoints(FunctionalBasis<Dim>* basis)
	{
		return CartesianShape<Dim>::GetNodalPoints(basis);
	}

	double SourceTerm(BasisFunction<Dim>* phi, SourceFunction* f)
	{
		function<double(RefPoint)> sourceTimesBasisFunction = [this, f, phi](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return f->Eval(domainPoint) * phi->Eval(refElementPoint);
		};

		return CartesianShape<Dim>::ComputeIntegral(sourceTimesBasisFunction, 0);
	}

	double L2ErrorPow2(function<double(RefPoint)> approximate, function<double(DomPoint)> exactSolution)
	{
		function<double(RefPoint)> errorFunction = [this, exactSolution, approximate](RefPoint refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return pow(exactSolution(domainPoint) - approximate(refElementPoint), 2);
		};

		return CartesianShape<Dim>::ComputeIntegral(errorFunction, 0);
	}

	//------------------------------------------------------------------//
	//                 Poisson_DG_Element implementation                //
	//------------------------------------------------------------------//
	
	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return CartesianShape<Dim>::MassTerm(phi1, phi2);
	}

	double StiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return CartesianShape<Dim>::StiffnessTerm(phi1, phi2);
	}

	//-------------------------------------------------------------------//
	//                 Poisson_HHO_Element implementation                //
	//-------------------------------------------------------------------//

	Eigen::MatrixXd CellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		return CartesianShape<Dim>::CellMassMatrix(basis);
	}

	Eigen::MatrixXd CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		return CartesianShape<Dim>::CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	double St(BasisFunction<Dim>* reconstructPhi1, BasisFunction<Dim>* reconstructPhi2)
	{
		return CartesianShape<Dim>::IntegralGradGradReconstruct(reconstructPhi1, reconstructPhi2);
	}
};