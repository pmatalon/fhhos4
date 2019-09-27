#pragma once
#include <Eigen/Sparse>
#include "../FunctionalBasis/FunctionalBasis.h"

template <int Dim>
class Poisson_HHO_Element;

template <int Dim>
class Poisson_HHO_Face : virtual public Face<Dim>
{
private:
	Eigen::MatrixXd _faceMassMatrix;
	Eigen::MatrixXd _invFaceMassMatrix;

	Eigen::MatrixXd _elem1_massCellFace;
	Eigen::MatrixXd _elem1_massReconstructFace;
	Eigen::MatrixXd _elem2_massCellFace;
	Eigen::MatrixXd _elem2_massReconstructFace;

	// Project a (k+1)-polynomial on the face
	Eigen::MatrixXd _elem1_projFromReconstruct;
	Eigen::MatrixXd _elem2_projFromReconstruct;

	// Project a k-polynomial on the face
	Eigen::MatrixXd _elem1_projFromCell;
	Eigen::MatrixXd _elem2_projFromCell;
public:
	HHOParameters<Dim>* HHO;

	Poisson_HHO_Face(BigNumber number, Element<Dim>* element1, Element<Dim>* element2) : Face<Dim>(number, element1, element2) {}

	void InitHHO(HHOParameters<Dim>* hho)
	{
		if (this->_faceMassMatrix.rows() > 0)
			return;

		this->HHO = hho;

		this->_faceMassMatrix = this->FaceMassMatrix(HHO->FaceBasis);
		this->_invFaceMassMatrix = this->_faceMassMatrix.inverse();

		this->_elem1_massReconstructFace = this->MassMatrix(HHO->FaceBasis, this->Element1, HHO->ReconstructionBasis);
		this->_elem1_massCellFace = this->MassMatrix(HHO->FaceBasis, this->Element1, HHO->CellBasis);
		this->_elem1_projFromReconstruct = this->_invFaceMassMatrix * _elem1_massReconstructFace;
		this->_elem1_projFromCell = this->_invFaceMassMatrix * _elem1_massCellFace;

		if (this->Element2 != NULL)
		{
			this->_elem2_massReconstructFace = this->MassMatrix(HHO->FaceBasis, this->Element2, HHO->ReconstructionBasis);
			this->_elem2_massCellFace = this->MassMatrix(HHO->FaceBasis, this->Element2, HHO->CellBasis);
			this->_elem2_projFromReconstruct = this->_invFaceMassMatrix * _elem2_massReconstructFace;
			this->_elem2_projFromCell = this->_invFaceMassMatrix * _elem2_massCellFace;
		}
	}

	Eigen::MatrixXd FaceMassMatrix()
	{
		return this->_faceMassMatrix;
	}

	Eigen::MatrixXd InvFaceMassMatrix()
	{
		return this->_invFaceMassMatrix;
	}

	virtual Eigen::MatrixXd FaceMassMatrix(FunctionalBasis<Dim-1>* basis) = 0;

	Eigen::MatrixXd MassMatrix(FunctionalBasis<Dim - 1>* basis, Element<Dim>* element, FunctionalBasis<Dim>* cellBasis)
	{
		Eigen::MatrixXd M(basis->LocalFunctions.size(), cellBasis->LocalFunctions.size());
		for (BasisFunction<Dim - 1>* phi1 : basis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : cellBasis->LocalFunctions)
			{
				double term = this->ComputeMassTerm(phi1, element, phi2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
			}
		}
		return M;
	}

	double ComputeMassTerm(BasisFunction<Dim - 1>* facePhi, Element<Dim>* element, BasisFunction<Dim>* reconstructPhi)
	{
		auto reconstructPhiOnFace = element->EvalPhiOnFace(this, reconstructPhi);

		function<double(RefPoint)> functionToIntegrate = [facePhi, reconstructPhiOnFace](RefPoint p) {
			return facePhi->Eval(p) * reconstructPhiOnFace(p);
		};

		int polynomialDegree = facePhi->GetDegree() + reconstructPhi->GetDegree();
		return this->ComputeIntegral(functionToIntegrate, polynomialDegree);
	}

	Eigen::MatrixXd GetMassCellFace(Element<Dim>* element)
	{
		if (element == this->Element1)
			return _elem1_massCellFace;
		else if (element == this->Element2)
			return _elem2_massCellFace;
		assert(false);
	}

	Eigen::MatrixXd GetMassReconstructFace(Element<Dim>* element)
	{
		if (element == this->Element1)
			return _elem1_massReconstructFace;
		else if (element == this->Element2)
			return _elem2_massReconstructFace;
		assert(false);
	}

	Eigen::MatrixXd GetProjFromReconstruct(Element<Dim>* element)
	{
		if (element == this->Element1)
			return _elem1_projFromReconstruct;
		else if (element == this->Element2)
			return _elem2_projFromReconstruct;
		assert(false);
	}

	Eigen::MatrixXd GetProjFromCell(Element<Dim>* element)
	{
		if (element == this->Element1)
			return _elem1_projFromCell;
		else if (element == this->Element2)
			return _elem2_projFromCell;
		assert(false);
	}

	Eigen::MatrixXd GetProjFromCell(Element<Dim>* element, FunctionalBasis<Dim>* cellInterpolationBasis)
	{
		Eigen::MatrixXd massFaceCell = this->MassMatrix(HHO->FaceBasis, element, cellInterpolationBasis);
		Eigen::MatrixXd projFromCell = this->_invFaceMassMatrix * massFaceCell;
		return projFromCell;
	}

	Eigen::MatrixXd ComputeCanonicalInjectionMatrixCoarseToFine(FunctionalBasis<Dim-1>* faceBasis)
	{
		Eigen::MatrixXd J(faceBasis->Size() * this->FinerFaces.size(), faceBasis->Size());

		for (auto f : this->FinerFaces)
		{
			Poisson_HHO_Face<Dim>* fineFace = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

			Eigen::MatrixXd fineCoarseMass(faceBasis->Size(), faceBasis->Size());
			for (BasisFunction<Dim-1>* finePhi : faceBasis->LocalFunctions)
			{
				for (BasisFunction<Dim-1>* coarsePhi : faceBasis->LocalFunctions)
				{
					function<double(RefPoint)> functionToIntegrate = [this, fineFace, finePhi, coarsePhi](RefPoint fineRefPoint) {
						DomPoint domPoint = fineFace->ConvertToDomain(fineRefPoint);
						RefPoint coarseRefPoint = this->ConvertToReference(domPoint);
						return finePhi->Eval(fineRefPoint)*coarsePhi->Eval(coarseRefPoint);
					};

					int polynomialDegree = finePhi->GetDegree() + coarsePhi->GetDegree();
					double integral = fineFace->ComputeIntegral(functionToIntegrate, polynomialDegree);
					fineCoarseMass(finePhi->LocalNumber, coarsePhi->LocalNumber) = integral;
				}
			}

			Eigen::MatrixXd invFineMass = fineFace->InvFaceMassMatrix();

			J.block(this->LocalNumberOf(fineFace)*faceBasis->Size(), 0, faceBasis->Size(), faceBasis->Size()) = invFineMass * fineCoarseMass;
		}

		return J;
	}
};