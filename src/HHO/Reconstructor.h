#pragma once
#include <Eigen/Sparse>
#include "Poisson_HHO_Element.h"
#include "Poisson_HHO_Face.h"

template <int Dim>
class Reconstructor
{
private:
	Poisson_HHO_Element<Dim>* _element;

	Eigen::MatrixXd _projFromReconstruct;
public:
	FunctionalBasis<Dim>* ReconstructionBasis;
	FunctionalBasis<Dim>* CellBasis;
	FunctionalBasis<Dim - 1>* FaceBasis;

	// Reconstruction operator as a matrix
	Eigen::MatrixXd P;

	// Consistency contribution
	Eigen::MatrixXd Acons;

	// Stabilization contribution
	Eigen::MatrixXd Astab;

	Reconstructor(Poisson_HHO_Element<Dim>* element, FunctionalBasis<Dim>* reconstructionBasis, FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim - 1>* faceBasis) :
		P(reconstructionBasis->Size(), cellBasis->Size() + element->Faces.size() * faceBasis->Size()),
		Acons(cellBasis->Size() + element->Faces.size() * faceBasis->Size(), cellBasis->Size() + element->Faces.size() * faceBasis->Size()),
		Astab(cellBasis->Size() + element->Faces.size() * faceBasis->Size(), cellBasis->Size() + element->Faces.size() * faceBasis->Size())
	{
		this->_element = element;
		this->ReconstructionBasis = reconstructionBasis;
		this->CellBasis = cellBasis;
		this->FaceBasis = faceBasis;

		Eigen::MatrixXd Mt = this->_element->MassMatrix(cellBasis);
		Eigen::MatrixXd Nt = this->_element->MassMatrix(cellBasis, reconstructionBasis);
		this->_projFromReconstruct = Mt.inverse() * Nt;

		for (auto f : this->_element->Faces)
		{
			Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
			face->InitHHO(reconstructionBasis, cellBasis, faceBasis);
		}

		this->AssembleReconstructionAndConsistencyMatrices();
		this->AssembleStabilizationMatrix();

		/*for (int i = 0; i < reconstructionBasis->Size(); i++)
		{
			Eigen::VectorXd vector(reconstructionBasis->Size());
			vector.setZero(reconstructionBasis->Size());
			vector(i) = 1;
			//cout << "------------- vector -------------" << endl << vector << endl;
			Eigen::VectorXd result = Reconstruct(Interpolate(vector));
			//cout << "------------- result -------------" << endl << result << endl;
		}*/
	}

	Eigen::MatrixXd ReconstructionMatrix()
	{
		return this->P;
	}

	Eigen::VectorXd Reconstruct(Eigen::VectorXd hybridVector)
	{
		return this->P * hybridVector;
	}

	Eigen::VectorXd Interpolate(Eigen::VectorXd vector)
	{
		Eigen::VectorXd hybridVector(this->CellBasis->Size() + this->_element->Faces.size() * this->FaceBasis->Size());
		Eigen::MatrixXd ProjT = _projFromReconstruct;

		hybridVector.head(this->CellBasis->Size()) = ProjT * vector;

		int index = this->CellBasis->Size();
		for (auto f : this->_element->Faces)
		{
			Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
			Eigen::MatrixXd ProjF = face->GetProjFromReconstruct(this->_element);

			hybridVector.segment(index, this->FaceBasis->Size()) = ProjF * vector;
			//cout << "------------- hybridVector -------------" << endl << hybridVector << endl;
			index += this->FaceBasis->Size();
		}
		return hybridVector;
	}

	double ConsistencyTerm(BasisFunction<Dim>* cellPhi1, BasisFunction<Dim>* cellPhi2)
	{
		return this->Acons(DOFNumber(cellPhi1), DOFNumber(cellPhi2));
	}
	double ConsistencyTerm(Face<Dim>* face, BasisFunction<Dim>* cellPhi, BasisFunction<Dim-1>* facePhi)
	{
		return this->Acons(DOFNumber(cellPhi), DOFNumber(face, facePhi));
	}
	double ConsistencyTerm(Face<Dim>* face1, BasisFunction<Dim - 1>* facePhi1, Face<Dim>* face2, BasisFunction<Dim - 1>* facePhi2)
	{
		return this->Acons(DOFNumber(face1, facePhi1), DOFNumber(face2, facePhi2));
	}

	double StabilizationTerm(BasisFunction<Dim>* cellPhi1, BasisFunction<Dim>* cellPhi2)
	{
		return this->Astab(DOFNumber(cellPhi1), DOFNumber(cellPhi2));
	}
	double StabilizationTerm(Face<Dim>* face, BasisFunction<Dim>* cellPhi, BasisFunction<Dim - 1>* facePhi)
	{
		return this->Astab(DOFNumber(cellPhi), DOFNumber(face, facePhi));
	}
	double StabilizationTerm(Face<Dim>* face1, BasisFunction<Dim - 1>* facePhi1, Face<Dim>* face2, BasisFunction<Dim - 1>* facePhi2)
	{
		return this->Astab(DOFNumber(face1, facePhi1), DOFNumber(face2, facePhi2));
	}

	double ReconstructionTerm(BasisFunction<Dim>* reconstructPhi, BasisFunction<Dim>* cellPhi2)
	{
		return this->P(reconstructPhi->LocalNumber, DOFNumber(cellPhi2));
	}
	double ReconstructionTerm(BasisFunction<Dim>* reconstructPhi, Face<Dim>* face, BasisFunction<Dim - 1>* facePhi)
	{
		return this->P(reconstructPhi->LocalNumber, DOFNumber(face, facePhi));
	}

private:

	//------------------------------------------------------------------------------------//
	//--------------- Reconstruction operator and consistency contribution ---------------//
	//------------------------------------------------------------------------------------//

	void AssembleReconstructionAndConsistencyMatrices()
	{
		Eigen::MatrixXd matrixToInvert = this->AssembleMatrixToInvert();
		Eigen::MatrixXd rhsMatrix = this->AssembleRHSMatrix();

		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver = matrixToInvert.colPivHouseholderQr();
		for (int j = 0; j < rhsMatrix.cols(); j++)
		{
			//cout << "----------------- col " << j << " ----------------" << endl;
			Eigen::VectorXd col = solver.solve(rhsMatrix.col(j));
			//cout << col << endl << endl;

			// We don't keep the last element of the column, which is the Lagrange multiplier
			//cout << "----------------- col " << j << " without Lagrange ----------------" << endl;
			auto colWithoutLagrangeMultiplier = col.head(this->ReconstructionBasis->Size());
			//cout << colWithoutLagrangeMultiplier << endl << endl;
			//cout << "----------------- P(" << j << ") ----------------" << endl;
			//cout << this->P.col(j) << endl << endl;
			this->P.col(j) = colWithoutLagrangeMultiplier;

			//cout << "----------------- P ----------------" << endl;
			//cout << this->P << endl << endl;
		}
		//cout << this->P << endl << endl;

		this->Acons = this->P.transpose() * matrixToInvert.topLeftCorner(this->ReconstructionBasis->Size(), this->ReconstructionBasis->Size()) * this->P;
	}

	Eigen::MatrixXd AssembleMatrixToInvert()
	{
		Eigen::MatrixXd matrixToInvert(this->ReconstructionBasis->Size() + 1, this->ReconstructionBasis->Size() + 1);
		//cout << matrixToInvert << endl << "----------------- matrixToInvert ----------------" << endl;
		this->AssembleSt(matrixToInvert);
		//cout << matrixToInvert << endl << "---------------- after AssembleSt -----------------" << endl;
		this->AssembleLt(matrixToInvert);
		//cout << matrixToInvert << endl << "---------------- after AssembleLt -----------------" << endl;
		matrixToInvert.bottomRightCorner<1, 1>() << 0;
		//cout << matrixToInvert << endl << "---------------- after 0 -----------------" << endl;
		return matrixToInvert;
	}

	void AssembleSt(Eigen::MatrixXd &matrixToInvert)
	{
		for (BasisFunction<Dim>* phi1 : this->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : this->ReconstructionBasis->LocalFunctions)
				matrixToInvert(phi1->LocalNumber, phi2->LocalNumber) = this->_element->St(phi1, phi2);
		}
	}

	void AssembleLt(Eigen::MatrixXd &matrixToInvert)
	{
		int last = this->ReconstructionBasis->NumberOfLocalFunctionsInElement(NULL);
		for (BasisFunction<Dim>* phi : this->ReconstructionBasis->LocalFunctions)
		{
			double Lt = this->_element->Lt(phi);
			matrixToInvert(last, phi->LocalNumber) = Lt;
			matrixToInvert(phi->LocalNumber, last) = Lt;
		}
	}

	void AssembleMt(Eigen::MatrixXd &rhsMatrix)
	{
		int last = this->ReconstructionBasis->NumberOfLocalFunctionsInElement(NULL);
		for (BasisFunction<Dim>* phi : this->CellBasis->LocalFunctions)
			rhsMatrix(last, phi->LocalNumber) = this->_element->Lt(phi);
	}

	Eigen::MatrixXd AssembleRHSMatrix()
	{
		auto nFaceUnknowns = this->_element->Faces.size() * this->FaceBasis->Size();
		auto nColumns = this->CellBasis->Size() + nFaceUnknowns;
		Eigen::MatrixXd rhsMatrix(this->ReconstructionBasis->Size() + 1, nColumns);
		//cout << rhsMatrix << endl << "----------------- rhsMatrix ----------------" << endl;
		
		// Top-left corner
		this->AssembleBt(rhsMatrix);
		//cout << rhsMatrix << endl << "----------------- after AssembleBt ----------------" << endl;
		
		// Top-right corner
		for (auto face : this->_element->Faces)
			this->AssembleBf(rhsMatrix, face);
		//cout << rhsMatrix << endl << "----------------- after AssembleBf ----------------" << endl;

		// Bottom-left corner
		this->AssembleMt(rhsMatrix);
		//cout << rhsMatrix << endl << "----------------- after AssembleMt ----------------" << endl;

		// Bottom-right corner
		rhsMatrix.bottomRightCorner(1, nFaceUnknowns) << Eigen::ArrayXXd::Zero(1, nFaceUnknowns);
		//cout << rhsMatrix << endl << "---------------------------------" << endl;

		return rhsMatrix;
	}

	void AssembleBt(Eigen::MatrixXd &rhsMatrix)
	{
		for (BasisFunction<Dim>* reconstructPhi : this->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim>* cellPhi : this->CellBasis->LocalFunctions)
				rhsMatrix(reconstructPhi->LocalNumber, DOFNumber(cellPhi)) = this->_element->Bt(reconstructPhi, cellPhi);
		}
	}

	void AssembleBf(Eigen::MatrixXd &rhsMatrix, Face<Dim>* face)
	{
		for (BasisFunction<Dim>* reconstructPhi : this->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim-1>* facePhi : this->FaceBasis->LocalFunctions)
				rhsMatrix(reconstructPhi->LocalNumber, DOFNumber(face, facePhi)) = this->_element->Bf(reconstructPhi, facePhi, face);
		}
	}

	//---------------------------------------------//
	//--------------- Stabilization ---------------//
	//---------------------------------------------//

	void AssembleStabilizationMatrix()
	{
		this->Astab = Eigen::MatrixXd::Zero(this->Astab.rows(), this->Astab.cols());

		Eigen::MatrixXd ProjT = _projFromReconstruct;
		//cout << "------------- ProjT -------------" << endl << ProjT << endl;
		Eigen::MatrixXd Dt = ProjT * this->P;
		//cout << "------------- Dt -------------" << endl << Dt << endl;
		for (int i = 0; i < Dt.rows(); i++)
			Dt(i, i) -= 1;
		//cout << "------------- Dt -------------" << endl << Dt << endl;

		for (auto f : this->_element->Faces)
		{
			Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
			Eigen::MatrixXd Mf = face->GetMassMatrix();
			Eigen::MatrixXd ProjF = face->GetProjFromReconstruct(this->_element);
			Eigen::MatrixXd ProjFT = face->GetProjFromCell(this->_element);

			Eigen::MatrixXd Df = ProjF * this->P;
			for (int i = 0; i < Df.rows(); i++)
				Df(i, FirstDOFNumber(face) + i) -= 1;

			Eigen::MatrixXd DiffTF = Df - ProjFT * Dt;
			double h = face->GetDiameter();
			this->Astab += DiffTF.transpose() * Mf * DiffTF / h;
			//cout << "------------- Astab -------------" << endl << Astab << endl;
		}
	}
	int DOFNumber(BasisFunction<Dim>* cellPhi)
	{
		return cellPhi->LocalNumber;
	}
	int DOFNumber(Face<Dim>* face, BasisFunction<Dim - 1>* facePhi)
	{
		return FirstDOFNumber(face) + facePhi->LocalNumber;
	}
public:
	int FirstDOFNumber(Face<Dim>* face)
	{
		return this->CellBasis->Size() + this->_element->LocalNumberOf(face) * this->FaceBasis->Size();
	}
};