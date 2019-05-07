#pragma once
#include <Eigen/Sparse>
#include "../FunctionalBasis/FunctionalBasis.h"

template <int Dim>
class Poisson_HHO_Element;

template <int Dim>
class Poisson_HHO_Face : virtual public Face<Dim>
{
private:
	Eigen::MatrixXd _massMatrix;

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
	Poisson_HHO_Face(BigNumber number, Element<Dim>* element1, Element<Dim>* element2) : Face<Dim>(number, element1, element2) {}

	void InitHHO(FunctionalBasis<Dim>* reconstructionBasis, FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim - 1>* faceBasis)
	{
		if (this->_massMatrix.rows() > 0)
			return;

		this->_massMatrix = this->MassMatrix(faceBasis);
		Eigen::MatrixXd invMf = this->_massMatrix.inverse();

		//cout << "------------- Mf -------------" << endl << Mf << endl;
		this->_elem1_massReconstructFace = this->MassMatrix(faceBasis, this->Element1, reconstructionBasis);
		//cout << "------------- Nf -------------" << endl << Nf << endl;
		this->_elem1_massCellFace = this->MassMatrix(faceBasis, this->Element1, cellBasis);
		//cout << "------------- Nft -------------" << endl << Nft << endl;
		this->_elem1_projFromReconstruct = invMf * _elem1_massReconstructFace;
		this->_elem1_projFromCell = invMf * _elem1_massCellFace;

		if (this->Element2 != NULL)
		{
			this->_elem2_massReconstructFace = this->MassMatrix(faceBasis, this->Element2, reconstructionBasis);
			//cout << "------------- Nf -------------" << endl << Nf << endl;
			this->_elem2_massCellFace = this->MassMatrix(faceBasis, this->Element2, cellBasis);
			//cout << "------------- Nft -------------" << endl << Nft << endl;
			this->_elem2_projFromReconstruct = invMf * _elem2_massReconstructFace;
			this->_elem2_projFromCell = invMf * _elem2_massCellFace;
		}
	}

	Eigen::MatrixXd GetMassMatrix()
	{
		return this->_massMatrix;
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

	Eigen::MatrixXd GetProjFromCoarserReconstruct()
	{
		Poisson_HHO_Element<Dim>* fineElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(this->Element1);

	}
};