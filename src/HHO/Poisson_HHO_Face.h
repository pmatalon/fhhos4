#pragma once
#include <Eigen/Sparse>
#include "../FunctionalBasis/FunctionalBasis.h"

template <int Dim>
class Poisson_HHO_Face : virtual public Face<Dim>
{
private:
	Eigen::MatrixXd _massMatrix;

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
		Eigen::MatrixXd Nf1 = this->MassMatrix(faceBasis, this->Element1, reconstructionBasis);
		//cout << "------------- Nf -------------" << endl << Nf << endl;
		Eigen::MatrixXd Nft1 = this->MassMatrix(faceBasis, this->Element1, cellBasis);
		//cout << "------------- Nft -------------" << endl << Nft << endl;
		this->_elem1_projFromReconstruct = invMf * Nf1;
		this->_elem1_projFromCell = invMf * Nft1;

		if (this->Element2 != NULL)
		{
			Eigen::MatrixXd Nf2 = this->MassMatrix(faceBasis, this->Element2, reconstructionBasis);
			//cout << "------------- Nf -------------" << endl << Nf << endl;
			Eigen::MatrixXd Nft2 = this->MassMatrix(faceBasis, this->Element2, cellBasis);
			//cout << "------------- Nft -------------" << endl << Nft << endl;
			this->_elem2_projFromReconstruct = invMf * Nf2;
			this->_elem2_projFromCell = invMf * Nft2;
		}
	}

	Eigen::MatrixXd GetMassMatrix()
	{
		return this->_massMatrix;
	}

	Eigen::MatrixXd GetProjFromReconstruct(Element<Dim>* element)
	{
		if (element == this->Element1)
			return _elem1_projFromReconstruct;
		else
			return _elem2_projFromReconstruct;
	}

	Eigen::MatrixXd GetProjFromCell(Element<Dim>* element)
	{
		if (element == this->Element1)
			return _elem1_projFromCell;
		else
			return _elem2_projFromCell;
	}

	Eigen::MatrixXd GetProjFromCoarserReconstruct()
	{
		Poisson_HHO_Element<Dim>* fineElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(this->Element1);

	}
};