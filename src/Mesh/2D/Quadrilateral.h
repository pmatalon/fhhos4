#pragma once
#include "../../DG/Poisson_DG_Element.h"
#include "../../HHO/Poisson_HHO_Element.h"
#include "QuadrilateralShape.h"
#include "Edge.h"
using namespace std;

class Quadrilateral : public Poisson_DG_Element<2>, public Poisson_HHO_Element<2>
{
private:
	QuadrilateralShape _shape;

public:
	Quadrilateral(int number, Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4) :
		Element(number),
		Poisson_DG_Element<2>(number),
		Poisson_HHO_Element<2>(number),
		_shape(v1, v2, v3, v4)
	{}

	inline Vertex* V1()
	{
		return _shape.V1;
	}
	inline Vertex* V2()
	{
		return _shape.V2;
	}
	inline Vertex* V3()
	{
		return _shape.V3;
	}
	inline Vertex* V4()
	{
		return _shape.V4;
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	const GeometricShapeWithReferenceShape<2>* Shape() const
	{
		return &_shape;
	}

	DimVector<2> OuterNormalVector(Face<2>* face)
	{
		DimVector<2> n;
		Edge* edge = dynamic_cast<Edge*>(face);
		Vertex* A = edge->Vertex1();
		Vertex* B = edge->Vertex2();

		// Condition 1: n.AB = 0
		// =>  n = (-AB.Y, AB.X)
		n << A->Y - B->Y, B->X - A->X;

		// Condition 2: n.AC < 0
		Vertex* C = nullptr;
		/*if (edge->Vertex1() == _shape.V1)
		{
			if (edge->Vertex2() == _shape.V2)
				C = _shape.V3;
			else
				C = _shape.V2;
		}
		else if (edge->Vertex1() == _shape.V2)
		{
			if (edge->Vertex2() == _shape.V1)
				C = _shape.V3;
			else
				C = _shape.V1;
		}
		else if (edge->Vertex1() == _shape.V3)
		{
			if (edge->Vertex2() == _shape.V1)
				C = _shape.V2;
			else
				C = _shape.V1;
		}
		else
			assert(false);*/

		DimVector<2> AC = Vect(A, C);
		//double nAC = n(0) * (C->X - A->X) + n(1) * (C->Y - A->Y);
		//if (nAC > 0)
		if (n.dot(AC) > 0)
			n = -1 * n;
		n = n.normalized();
		return n;
	}

	//-------------------------------------------------------------------//
	//                  Poisson_DG_Element implementation                //
	//-------------------------------------------------------------------//

	double MassTerm(BasisFunction<2>* phi1, BasisFunction<2>* phi2)
	{
		return _shape.MassTerm(phi1, phi2);
	}

	double StiffnessTerm(BasisFunction<2>* phi1, BasisFunction<2>* phi2)
	{
		Tensor<2>* K = Tensor<2>::Isotropic();
		double integral = Shape()->ComputeIntegralKGradGrad(K, phi1, phi2);
		delete K;
		return integral;
	}

	//-------------------------------------------------------------------//
	//                 Poisson_HHO_Element implementation                //
	//-------------------------------------------------------------------//

	DenseMatrix CellMassMatrix(FunctionalBasis<2>* basis)
	{
		return _shape.CellMassMatrix(basis);
	}

	DenseMatrix CellReconstructMassMatrix(FunctionalBasis<2>* cellBasis, FunctionalBasis<2>* reconstructBasis)
	{
		return _shape.CellReconstructMassMatrix(cellBasis, reconstructBasis);
	}

	double IntegralKGradGradReconstruct(Tensor<2>* K, BasisFunction<2>* reconstructPhi1, BasisFunction<2>* reconstructPhi2)
	{
		return Shape()->ComputeIntegralKGradGrad(K, reconstructPhi1, reconstructPhi2);
	}

};