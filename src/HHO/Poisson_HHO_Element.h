#pragma once
#include "../Utils/SourceFunction.h"
#include "../FunctionalBasis/FunctionalBasis.h"
#include "../Mesh/Face.h"

template <short Dim>
class Poisson_HHO_Element : virtual public Element<Dim>
{
public:
	Poisson_HHO_Element(BigNumber number) : Element<Dim>(number) {}

	virtual double St(BasisFunction<Dim>* reconstructPhi1, BasisFunction<Dim>* reconstructPhi2) = 0;
	virtual double Lt(BasisFunction<Dim>* phi) = 0;
	virtual double Bt(BasisFunction<Dim>* reconstructPhi, BasisFunction<Dim>* cellPhi) = 0;
	virtual double Bf(BasisFunction<Dim>* reconstructPhi, BasisFunction<Dim-1>* facePhi, Face<Dim>* face) = 0;
	
	virtual void InitReconstructor(FunctionalBasis<Dim>* reconstructionBasis, FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim - 1>* faceBasis) = 0;
	virtual Eigen::VectorXd Reconstruct(Eigen::VectorXd hybridVector) = 0;

	virtual int FirstDOFLocalNumber(Face<Dim>* face) = 0;

	virtual double ConsistencyTerm(BasisFunction<Dim>* cellPhi1, BasisFunction<Dim>* cellPhi2) = 0;
	virtual double ConsistencyTerm(Face<Dim>* face, BasisFunction<Dim>* cellPhi, BasisFunction<Dim-1>* facePhi) = 0;
	virtual double ConsistencyTerm(Face<Dim>* face1, BasisFunction<Dim-1>* facePhi1, Face<Dim>* face2, BasisFunction<Dim - 1>* facePhi2) = 0;

	virtual double StabilizationTerm(BasisFunction<Dim>* cellPhi1, BasisFunction<Dim>* cellPhi2) = 0;
	virtual double StabilizationTerm(Face<Dim>* face, BasisFunction<Dim>* cellPhi, BasisFunction<Dim - 1>* facePhi) = 0;
	virtual double StabilizationTerm(Face<Dim>* face1, BasisFunction<Dim - 1>* facePhi1, Face<Dim>* face2, BasisFunction<Dim - 1>* facePhi2) = 0;

	virtual double ReconstructionTerm(BasisFunction<Dim>* reconstrucPhi, BasisFunction<Dim>* cellPhi) = 0;
	virtual double ReconstructionTerm(BasisFunction<Dim>* reconstrucPhi, Face<Dim>* face, BasisFunction<Dim-1>* facePhi) = 0;

	virtual double SourceTerm(BasisFunction<Dim>* cellPhi, SourceFunction* f) = 0;


	virtual ~Poisson_HHO_Element() {}
};