#pragma once
#include "../../DG/Poisson_DG_Element.h"
#include "../../HHO/Poisson_HHO_Element.h"
#include "GenericEdge.h"
#include "../ReferenceTriangle.h"
using namespace std;

class Triangle : public Poisson_DG_Element<2>, public Poisson_HHO_Element<2>
{
private:
	double _diameter;
	double _measure;
	DomPoint _center;

	DimMatrix<2> _inverseJacobianTranspose;
	double _detJacobian;

public:
	Vertex* V1;
	Vertex* V2;
	Vertex* V3;

	static ReferenceTriangle RefTriangle;

	Triangle(int number, Vertex* v1, Vertex* v2, Vertex* v3) :
		Element(number),
		Poisson_DG_Element<2>(number),
		Poisson_HHO_Element<2>(number)
	{
		V1 = v1;
		V2 = v2;
		V3 = v3;
		Init();
	}

	void Init()
	{
		double lengthEdge12 = sqrt(pow(V2->X - V1->X, 2) + pow(V2->Y - V1->Y, 2));
		double lengthEdge23 = sqrt(pow(V3->X - V2->X, 2) + pow(V3->Y - V2->Y, 2));
		double lengthEdge13 = sqrt(pow(V3->X - V1->X, 2) + pow(V3->Y - V1->Y, 2));
		_diameter = max(lengthEdge12, max(lengthEdge23, lengthEdge13));

		_measure = 0.5 * abs(V1->X * (V2->Y - V3->Y) + V2->X * (V3->Y - V1->Y) + V3->X * (V1->Y - V2->Y));

		_center = DomPoint((V1->X + V2->X + V3->X) / 3, (V1->Y + V2->Y + V3->Y) / 3);

		_detJacobian = _measure / RefTriangle.Measure();

		DimMatrix<2> inverseJacobian;
		inverseJacobian(0, 0) =  (V3->Y - V1->Y) / ((V3->Y - V1->Y)*(V2->X - V1->X) - (V3->X - V1->X)*(V2->Y - V1->Y));
		inverseJacobian(0, 1) = -(V3->X - V1->X) / ((V3->Y - V1->Y)*(V2->X - V1->X) - (V3->X - V1->X)*(V2->Y - V1->Y));
		inverseJacobian(1, 0) =  (V2->Y - V1->Y) / ((V2->Y - V1->Y)*(V3->X - V1->X) - (V2->X - V1->X)*(V3->Y - V1->Y));
		inverseJacobian(1, 1) = -(V2->X - V1->X) / ((V2->Y - V1->Y)*(V3->X - V1->X) - (V2->X - V1->X)*(V3->Y - V1->Y));
		_inverseJacobianTranspose = inverseJacobian.transpose();
	}

	void Serialize(ostream& os) const override
	{
		Element<2>::Serialize(os);
		os << ", Triangle";
		os << " ";
		V1->Serialize(os, 2);
		os << "--";
		V2->Serialize(os, 2);
		os << "--";
		V3->Serialize(os, 2);
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	inline double Diameter()
	{
		return _diameter;
	}
	inline double Measure()
	{
		return _measure;
	}
	inline DomPoint Center()
	{
		return _center;
	}

	DimVector<2> OuterNormalVector(Face<2>* face)
	{
		DimVector<2> n;
		Edge* edge = dynamic_cast<Edge*>(face);
		Vertex* A = edge->Vertex1;
		Vertex* B = edge->Vertex2;

		// Condition 1: n.AB = 0
		// =>  n = (-AB.Y, AB.X)
		n << A->Y - B->Y, B->X - A->X;

		// Condition 2: n.AC < 0
		Vertex* C = nullptr;
		if (edge->Vertex1 == this->V1)
		{
			if (edge->Vertex2 == this->V2)
				C = V3;
			else
				C = V2;
		}
		else if (edge->Vertex1 == this->V2)
		{
			if (edge->Vertex2 == this->V1)
				C = V3;
			else
				C = V1;
		}
		else if (edge->Vertex1 == this->V3)
		{
			if (edge->Vertex2 == this->V1)
				C = V2;
			else
				C = V1;
		}
		else
			assert(false);

		DimVector<2> AC = Vect(A, C);
		//double nAC = n(0) * (C->X - A->X) + n(1) * (C->Y - A->Y);
		//if (nAC > 0)
		if (n.dot(AC) > 0)
			n = -1 * n;
		n = n.normalized();
		return n;
	}

	inline double DetJacobian() const
	{
		return _detJacobian;
	}
	inline DimMatrix<2> InverseJacobianTranspose() const
	{
		return _inverseJacobianTranspose;
	}

	DomPoint ConvertToDomain(RefPoint refPoint) const
	{
		double t = refPoint.X;
		double u = refPoint.Y;

		DomPoint p;
		p.X = (V2->X - V1->X) * t + (V3->X - V1->X)*u + V1->X;
		p.Y = (V2->Y - V1->Y) * t + (V3->Y - V1->Y)*u + V1->Y;
		return p;
	}

	RefPoint ConvertToReference(DomPoint domainPoint) const
	{
		double x = domainPoint.X;
		double y = domainPoint.Y;

		double t = ((V3->Y - V1->Y)*(x - V1->X) - (V3->X - V1->X)*(y - V1->Y)) / ((V3->Y - V1->Y)*(V2->X - V1->X) - (V3->X - V1->X)*(V2->Y - V1->Y));
		double u = ((V2->Y - V1->Y)*(x - V1->X) - (V2->X - V1->X)*(y - V1->Y)) / ((V2->Y - V1->Y)*(V3->X - V1->X) - (V2->X - V1->X)*(V3->Y - V1->Y));
		RefPoint p(t, u);
		return p;
	}

	double Integral(BasisFunction<2>* phi) const
	{
		return DetJacobian() * RefTriangle.Integral(phi);
	}
	double Integral(RefFunction func) const
	{
		return DetJacobian() * RefTriangle.Integral(func);
	}
	double Integral(RefFunction func, int polynomialDegree) const
	{
		return DetJacobian() * RefTriangle.Integral(func, polynomialDegree);
	}

	//-------------------------------------------------------------------//
	//                  Poisson_DG_Element implementation                //
	//-------------------------------------------------------------------//

	double MassTerm(BasisFunction<2>* phi1, BasisFunction<2>* phi2)
	{
		return DetJacobian() * RefTriangle.MassTerm(phi1, phi2);
	}

	double StiffnessTerm(BasisFunction<2>* phi1, BasisFunction<2>* phi2)
	{
		Tensor<2>* K = Tensor<2>::Isotropic();
		double integral = ComputeIntegralKGradGrad(K, phi1, phi2);
		delete K;
		return integral;
	}

	//-------------------------------------------------------------------//
	//                 Poisson_HHO_Element implementation                //
	//-------------------------------------------------------------------//

	DenseMatrix CellMassMatrix(FunctionalBasis<2>* basis)
	{
		return DetJacobian() * RefTriangle.CellMassMatrix(basis);
	}

	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<2>* cellBasis, FunctionalBasis<2>* reconstructBasis)
	{
		return DetJacobian() * RefTriangle.CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	double IntegralKGradGradReconstruct(Tensor<2>* K, BasisFunction<2>* reconstructPhi1, BasisFunction<2>* reconstructPhi2)
	{
		return ComputeIntegralKGradGrad(K, reconstructPhi1, reconstructPhi2);
	}

	double ComputeIntegralKGradGrad(Tensor<2>* K, BasisFunction<2>* phi1, BasisFunction<2>* phi2) const
	{
		if (phi1->GetDegree() == 0 || phi2->GetDegree() == 0)
			return 0;

		DimMatrix<2> invJ = InverseJacobianTranspose();

		RefFunction functionToIntegrate = [K, phi1, phi2, invJ](RefPoint p) {
			DimVector<2> gradPhi1 = invJ * phi1->Grad(p);
			DimVector<2> gradPhi2 = invJ * phi2->Grad(p);
			return (K * gradPhi1).dot(gradPhi2);
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Integral(functionToIntegrate, polynomialDegree);
	}

};

ReferenceTriangle Triangle::RefTriangle = ReferenceTriangle();