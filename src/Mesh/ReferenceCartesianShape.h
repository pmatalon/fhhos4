#pragma once
#include <Eigen/Dense>
#include "../FunctionalBasis/FunctionalBasis.h"

template <int Dim>
class ReferenceCartesianShape
{
private:
	// For DG
	Eigen::MatrixXd _stiffnessMatrix;
	Eigen::MatrixXd _massMatrix;

	// For HHO
	Eigen::MatrixXd _cellMassMatrix;
	Eigen::MatrixXd _reconstructMassMatrix;
	Eigen::MatrixXd _cellStiffnessMatrix;
	Eigen::MatrixXd _reconstructStiffnessMatrix;
	Eigen::MatrixXd _faceMassMatrix;
	Eigen::MatrixXd _cellReconstructMassMatrix;
public:
	ReferenceCartesianShape()
	{
	}

	double MassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->_massMatrix(phi1->LocalNumber, phi2->LocalNumber);
	}

	double StiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->_stiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber);
	}

	double ReconstructStiffnessTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		return this->_reconstructStiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber);
	}

	double ComputeIntegral(BasisFunction<Dim>* phi)
	{
		return Utils::Integral(phi);
	}

	double ComputeIntegral(function<double(RefPoint)> func, int polynomialDegree)
	{
		return Utils::Integral<Dim>(func, polynomialDegree);
	}

	double ComputeIntegral(function<double(RefPoint)> func)
	{
		return Utils::Integral<Dim>(func);
	}

	void ComputeAndStoreMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_massMatrix.rows() == 0)
			_massMatrix = ComputeAndReturnMassMatrix(basis);
	}
	void ComputeAndStoreCellMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_cellMassMatrix.rows() == 0)
			_cellMassMatrix = ComputeAndReturnMassMatrix(basis);
	}
	void ComputeAndStoreFaceMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_faceMassMatrix.rows() == 0)
			_faceMassMatrix = ComputeAndReturnMassMatrix(basis);
	}
	void ComputeAndStoreReconstructMassMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_reconstructMassMatrix.rows() == 0)
			_reconstructMassMatrix = ComputeAndReturnMassMatrix(basis);
	}
	void ComputeAndStoreStiffnessMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_stiffnessMatrix.rows() == 0)
			_stiffnessMatrix = ComputeAndReturnStiffnessMatrix(basis);
	}
	void ComputeAndStoreCellStiffnessMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_cellStiffnessMatrix.rows() == 0)
			_cellStiffnessMatrix = ComputeAndReturnStiffnessMatrix(basis);
	}
	void ComputeAndStoreReconstructStiffnessMatrix(FunctionalBasis<Dim>* basis)
	{
		if (_reconstructStiffnessMatrix.rows() == 0)
			_reconstructStiffnessMatrix = ComputeAndReturnStiffnessMatrix(basis);
	}
	void ComputeAndStoreCellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis)
	{
		if (_cellReconstructMassMatrix.rows() == 0)
			_cellReconstructMassMatrix = ComputeAndReturnMassMatrix(cellBasis, reconstructBasis);
	}

	Eigen::MatrixXd StoredCellMassMatrix()
	{
		return _cellMassMatrix;
	}
	Eigen::MatrixXd StoredReconstructMassMatrix()
	{
		return _reconstructMassMatrix;
	}
	Eigen::MatrixXd StoredFaceMassMatrix()
	{
		return _faceMassMatrix;
	}
	Eigen::MatrixXd StoredCellReconstructMassMatrix()
	{
		return _cellReconstructMassMatrix;
	}

private:
	Eigen::MatrixXd ComputeAndReturnMassMatrix(FunctionalBasis<Dim>* basis)
	{
		Eigen::MatrixXd M = Eigen::MatrixXd(basis->Size(), basis->Size());
		for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
			{
				if (phi2->LocalNumber > phi1->LocalNumber)
					break;
				double term = ComputeMassTerm(phi1, phi2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
				M(phi2->LocalNumber, phi1->LocalNumber) = term;
			}
		}
		return M;
	}

	Eigen::MatrixXd ComputeAndReturnMassMatrix(FunctionalBasis<Dim>* basis1, FunctionalBasis<Dim>* basis2)
	{
		Eigen::MatrixXd M(basis1->LocalFunctions.size(), basis2->LocalFunctions.size());
		for (BasisFunction<Dim>* phi1 : basis1->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : basis2->LocalFunctions)
			{
				double term = ComputeMassTerm(phi1, phi2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
			}
		}
		return M;
	}

	Eigen::MatrixXd ComputeAndReturnStiffnessMatrix(FunctionalBasis<Dim>* basis)
	{
		Eigen::MatrixXd stiffnessMatrix = Eigen::MatrixXd(basis->Size(), basis->Size());
		for (BasisFunction<Dim>* phi1 : basis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : basis->LocalFunctions)
				stiffnessMatrix(phi1->LocalNumber, phi2->LocalNumber) = ComputeIntegralGradGrad(phi1, phi2);
		}
		return stiffnessMatrix;
	}
public:
	double ComputeMassTerm(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		function<double(Point)> functionToIntegrate = [phi1, phi2](RefPoint p) {
			return phi1->Eval(p)*phi2->Eval(p);
		};

		int polynomialDegree = phi1->GetDegree() + phi2->GetDegree();
		return Utils::Integral<Dim>(functionToIntegrate, polynomialDegree);
	}

	double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2)
	{
		if (phi1->GetDegree() == 0 || phi1->GetDegree() == 0)
			return 0;

		function<double(Point)> functionToIntegrate = [phi1, phi2](RefPoint p) {
			return Element<Dim>::InnerProduct(phi1->Grad(p), phi2->Grad(p));
		};

		int polynomialDegree = max(0, phi1->GetDegree() + phi2->GetDegree() - 2);
		return Utils::Integral<Dim>(functionToIntegrate, polynomialDegree);
	}
};