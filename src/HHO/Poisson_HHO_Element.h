#pragma once
#include <Eigen/Sparse>
#include "../Utils/SourceFunction.h"
#include "../FunctionalBasis/FunctionalBasis.h"
#include "Poisson_HHO_Face.h"
#include "../Mesh/Face.h"

template <int Dim>
class Poisson_HHO_Element : virtual public Element<Dim>
{
private:
	Eigen::MatrixXd _cellMassMatrix;
	Eigen::MatrixXd _projFromReconstruct;
public:
	FunctionalBasis<Dim>* ReconstructionBasis;
	FunctionalBasis<Dim>* CellBasis;
	FunctionalBasis<Dim - 1>* FaceBasis;
	double Kappa; // constant diffusion coefficient

	// Reconstruction operator as a matrix
	Eigen::MatrixXd P;

	// Consistency contribution
	Eigen::MatrixXd Acons;

	// Stabilization contribution
	Eigen::MatrixXd Astab;

	// Local operator matrix (= Acons + Astab)
	Eigen::MatrixXd A;

	Eigen::MatrixXd invAtt;

	Poisson_HHO_Element(BigNumber number) : Element<Dim>(number) {}
	
	void InitHHO(FunctionalBasis<Dim>* reconstructionBasis, FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim - 1> * faceBasis, DiffusionPartition diffusionPartition)
	{
		this->ReconstructionBasis = reconstructionBasis;
		this->CellBasis = cellBasis;
		this->FaceBasis = faceBasis;
		this->Kappa = this->DiffusionCoefficient(diffusionPartition);

		this->_cellMassMatrix = this->CellMassMatrix(cellBasis);
		Eigen::MatrixXd Nt = this->CellReconstructMassMatrix(cellBasis, reconstructionBasis);
		this->_projFromReconstruct = _cellMassMatrix.inverse() * Nt;

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

		int nCellUnknowns = cellBasis->Size();
		int nTotalFaceUnknowns = this->Faces.size() * faceBasis->Size();

		this->A = Acons + Astab;
		auto Att = A.topLeftCorner(nCellUnknowns, nCellUnknowns);
		//auto Aff = A.bottomRightCorner(nTotalFaceUnknowns, nTotalFaceUnknowns);
		auto Atf = A.topRightCorner(nCellUnknowns, nTotalFaceUnknowns);
		this->invAtt = Att.inverse();
		
	}

	Eigen::MatrixXd CellMassMatrix()
	{
		return this->_cellMassMatrix;
	}
	
	Eigen::MatrixXd ComputeCanonicalInjectionMatrixCoarseToFine(FunctionalBasis<Dim>* cellBasis)
	{
		Eigen::MatrixXd J(cellBasis->Size() * this->FinerElements.size(), cellBasis->Size());

		for (auto e : this->FinerElements)
		{
			Poisson_HHO_Element<Dim>* fineElement = dynamic_cast<Poisson_HHO_Element<Dim>*>(e);

			Eigen::MatrixXd fineCoarseMass(cellBasis->Size(), cellBasis->Size());
			for (BasisFunction<Dim>* finePhi : cellBasis->LocalFunctions)
			{
				for (BasisFunction<Dim>* coarsePhi : cellBasis->LocalFunctions)
				{
					function<double(RefPoint)> functionToIntegrate = [this, fineElement, finePhi, coarsePhi](RefPoint fineRefPoint) {
						DomPoint domPoint = fineElement->ConvertToDomain(fineRefPoint);
						RefPoint coarseRefPoint = this->ConvertToReference(domPoint);
						return finePhi->Eval(fineRefPoint)*coarsePhi->Eval(coarseRefPoint);
					};
					
					int polynomialDegree = finePhi->GetDegree() + coarsePhi->GetDegree();
					double integral = fineElement->ComputeIntegral(functionToIntegrate, polynomialDegree);
					fineCoarseMass(finePhi->LocalNumber, coarsePhi->LocalNumber) = integral;
				}
			}
			
			Eigen::MatrixXd fineMass = fineElement->CellMassMatrix(cellBasis);

			J.block(this->LocalNumberOf(fineElement)*cellBasis->Size(), 0, cellBasis->Size(), cellBasis->Size()) = fineMass.inverse() * fineCoarseMass;
		}

		return J;
	}

	Eigen::VectorXd Reconstruct(Eigen::VectorXd hybridVector)
	{
		return this->P * hybridVector;
	}
	
	Eigen::MatrixXd ReconstructionMatrix()
	{
		return this->P;
	}

	Eigen::MatrixXd ReconstructionFromFacesMatrix()
	{
		int nCellUnknowns = CellBasis->Size();
		int nTotalFaceUnknowns = this->Faces.size() * FaceBasis->Size();

		auto Atf = this->A.topRightCorner(nCellUnknowns, nTotalFaceUnknowns);
		Eigen::MatrixXd solveCellUnknowns = -this->invAtt * Atf;

		Eigen::MatrixXd createHybridVectorFromFacesMatrix(nCellUnknowns + nTotalFaceUnknowns, nTotalFaceUnknowns);
		createHybridVectorFromFacesMatrix.topRows(nCellUnknowns) = solveCellUnknowns;
		createHybridVectorFromFacesMatrix.bottomRows(nTotalFaceUnknowns) = Eigen::MatrixXd::Identity(nTotalFaceUnknowns, nTotalFaceUnknowns);
		return this->P * createHybridVectorFromFacesMatrix;
	}

	/*Eigen::VectorXd Interpolate(Eigen::VectorXd reconstructVector)
	{
		Eigen::VectorXd hybridVector(this->CellBasis->Size() + this->Faces.size() * this->FaceBasis->Size());
		Eigen::MatrixXd ProjT = _projFromReconstruct;

		hybridVector.head(this->CellBasis->Size()) = ProjT * reconstructVector;

		int index = this->CellBasis->Size();
		for (auto f : this->Faces)
		{
			Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
			Eigen::MatrixXd ProjF = face->GetProjFromReconstruct(this);

			hybridVector.segment(index, this->FaceBasis->Size()) = ProjF * reconstructVector;
			index += this->FaceBasis->Size();
		}
		return hybridVector;
	}*/

	inline double MatrixTerm(BasisFunction<Dim>* cellPhi1, BasisFunction<Dim>* cellPhi2)
	{
		return this->A(DOFNumber(cellPhi1), DOFNumber(cellPhi2));
	}
	inline double MatrixTerm(Face<Dim>* face, BasisFunction<Dim>* cellPhi, BasisFunction<Dim - 1> * facePhi)
	{
		return this->A(DOFNumber(cellPhi), DOFNumber(face, facePhi));
	}
	inline double MatrixTerm(Face<Dim>* face1, BasisFunction<Dim - 1> * facePhi1, Face<Dim>* face2, BasisFunction<Dim - 1> * facePhi2)
	{
		return this->A(DOFNumber(face1, facePhi1), DOFNumber(face2, facePhi2));
	}

	inline double ConsistencyTerm(BasisFunction<Dim>* cellPhi1, BasisFunction<Dim>* cellPhi2)
	{
		return this->Acons(DOFNumber(cellPhi1), DOFNumber(cellPhi2));
	}
	inline double ConsistencyTerm(Face<Dim>* face, BasisFunction<Dim>* cellPhi, BasisFunction<Dim - 1> * facePhi)
	{
		return this->Acons(DOFNumber(cellPhi), DOFNumber(face, facePhi));
	}
	inline double ConsistencyTerm(Face<Dim>* face1, BasisFunction<Dim - 1> * facePhi1, Face<Dim>* face2, BasisFunction<Dim - 1> * facePhi2)
	{
		return this->Acons(DOFNumber(face1, facePhi1), DOFNumber(face2, facePhi2));
	}

	inline double StabilizationTerm(BasisFunction<Dim>* cellPhi1, BasisFunction<Dim>* cellPhi2)
	{
		return this->Astab(DOFNumber(cellPhi1), DOFNumber(cellPhi2));
	}
	inline double StabilizationTerm(Face<Dim>* face, BasisFunction<Dim>* cellPhi, BasisFunction<Dim - 1> * facePhi)
	{
		return this->Astab(DOFNumber(cellPhi), DOFNumber(face, facePhi));
	}
	inline double StabilizationTerm(Face<Dim>* face1, BasisFunction<Dim - 1> * facePhi1, Face<Dim>* face2, BasisFunction<Dim - 1> * facePhi2)
	{
		return this->Astab(DOFNumber(face1, facePhi1), DOFNumber(face2, facePhi2));
	}

	inline double ReconstructionTerm(BasisFunction<Dim>* reconstructPhi, BasisFunction<Dim>* cellPhi2)
	{
		return this->P(reconstructPhi->LocalNumber, DOFNumber(cellPhi2));
	}
	inline double ReconstructionTerm(BasisFunction<Dim>* reconstructPhi, Face<Dim>* face, BasisFunction<Dim - 1> * facePhi)
	{
		return this->P(reconstructPhi->LocalNumber, DOFNumber(face, facePhi));
	}


	virtual double SourceTerm(BasisFunction<Dim>* cellPhi, SourceFunction* f) = 0;
	virtual double IntegralGradGradReconstruct(BasisFunction<Dim>* reconstructPhi1, BasisFunction<Dim>* reconstructPhi2) = 0;
	virtual double ComputeIntegralGradGrad(BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) = 0;
	virtual Eigen::MatrixXd CellMassMatrix(FunctionalBasis<Dim>* basis) = 0;
	virtual Eigen::MatrixXd CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis) = 0;


	virtual ~Poisson_HHO_Element() {}

private:

	//------------------------------------------------------------------------------------//
	//--------------- Reconstruction operator and consistency contribution ---------------//
	//------------------------------------------------------------------------------------//

	void AssembleReconstructionAndConsistencyMatrices()
	{
		this->P = Eigen::MatrixXd(ReconstructionBasis->Size(), CellBasis->Size() + this->Faces.size() * FaceBasis->Size());

		Eigen::MatrixXd reconstructionMatrixToInvert = this->AssembleReconstructionMatrixToInvert();
		Eigen::MatrixXd rhsMatrix = this->AssembleRHSMatrix();

		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver = reconstructionMatrixToInvert.colPivHouseholderQr();
		for (int j = 0; j < rhsMatrix.cols(); j++)
		{
			Eigen::VectorXd col = solver.solve(rhsMatrix.col(j));

			// We don't keep the last element of the column, which is the Lagrange multiplier
			auto colWithoutLagrangeMultiplier = col.head(this->ReconstructionBasis->Size());
			this->P.col(j) = colWithoutLagrangeMultiplier;
		}

		Eigen::MatrixXd laplacianMatrix = reconstructionMatrixToInvert.topLeftCorner(this->ReconstructionBasis->Size(), this->ReconstructionBasis->Size()); // Grad-Grad part (block S)

		this->Acons = this->P.transpose() * laplacianMatrix * this->P;
	}

	Eigen::MatrixXd AssembleReconstructionMatrixToInvert()
	{
		Eigen::MatrixXd matrixToInvert(this->ReconstructionBasis->Size() + 1, this->ReconstructionBasis->Size() + 1);
		this->AssembleGradientReconstructionMatrix(matrixToInvert); // Block S
		this->AssembleMeanValueCondition(matrixToInvert); // Blocks L and L_transpose
		matrixToInvert.bottomRightCorner<1, 1>() << 0;
		return matrixToInvert;
	}

	void AssembleGradientReconstructionMatrix(Eigen::MatrixXd & reconstructionMatrixToInvert)
	{
		for (BasisFunction<Dim>* phi1 : this->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : this->ReconstructionBasis->LocalFunctions)
				reconstructionMatrixToInvert(phi1->LocalNumber, phi2->LocalNumber) = Kappa * this->IntegralGradGradReconstruct(phi1, phi2);
		}
	}

	void AssembleMeanValueCondition(Eigen::MatrixXd & reconstructionMatrixToInvert)
	{
		int last = this->ReconstructionBasis->NumberOfLocalFunctionsInElement(NULL);
		for (BasisFunction<Dim>* phi : this->ReconstructionBasis->LocalFunctions)
		{
			double meanValue = this->Integral(phi);
			reconstructionMatrixToInvert(last, phi->LocalNumber) = meanValue;
			reconstructionMatrixToInvert(phi->LocalNumber, last) = meanValue;
		}
	}

	void AssembleMeanValueConditionRHS(Eigen::MatrixXd & rhsMatrix)
	{
		int last = this->ReconstructionBasis->NumberOfLocalFunctionsInElement(NULL);
		for (BasisFunction<Dim>* phi : this->CellBasis->LocalFunctions)
			rhsMatrix(last, phi->LocalNumber) = this->Integral(phi);
	}

	Eigen::MatrixXd AssembleRHSMatrix()
	{
		auto nFaceUnknowns = this->Faces.size() * this->FaceBasis->Size();
		auto nColumns = this->CellBasis->Size() + nFaceUnknowns;
		Eigen::MatrixXd rhsMatrix(this->ReconstructionBasis->Size() + 1, nColumns);

		// Top-left corner (Block Bt)
		this->AssembleIntegrationByPartsRHS_cell(rhsMatrix);

		// Top-right corner (Block B_frontier)
		for (auto face : this->Faces)
			this->AssembleIntegrationByPartsRHS_face(rhsMatrix, face);

		// Bottom-left corner (Block Lt)
		this->AssembleMeanValueConditionRHS(rhsMatrix);

		// Bottom-right corner (0)
		rhsMatrix.bottomRightCorner(1, nFaceUnknowns) << Eigen::ArrayXXd::Zero(1, nFaceUnknowns);

		return rhsMatrix;
	}

	void AssembleIntegrationByPartsRHS_cell(Eigen::MatrixXd & rhsMatrix)
	{
		for (BasisFunction<Dim>* reconstructPhi : this->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim>* cellPhi : this->CellBasis->LocalFunctions)
				rhsMatrix(reconstructPhi->LocalNumber, DOFNumber(cellPhi)) = this->IntegrationByPartsRHS_cell(reconstructPhi, cellPhi);
		}
	}

	void AssembleIntegrationByPartsRHS_face(Eigen::MatrixXd & rhsMatrix, Face<Dim> * f)
	{
		Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
		for (BasisFunction<Dim>* reconstructPhi : this->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim - 1> * facePhi : this->FaceBasis->LocalFunctions)
				rhsMatrix(reconstructPhi->LocalNumber, DOFNumber(face, facePhi)) = this->IntegrationByPartsRHS_face(face, reconstructPhi, facePhi);
		}
	}

	double IntegrationByPartsRHS_cell(BasisFunction<Dim>* reconstructPhi, BasisFunction<Dim>* cellPhi)
	{
		if (reconstructPhi->GetDegree() == 0)
			return 0;

		double integralGradGrad = Kappa * this->ComputeIntegralGradGrad(reconstructPhi, cellPhi);

		double sumFaces = 0;
		for (auto f : this->Faces)
		{
			Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

			auto phi = this->EvalPhiOnFace(face, cellPhi);
			auto gradPhi = this->GradPhiOnFace(face, reconstructPhi);
			auto normal = this->OuterNormalVector(face);

			std::function<double(RefPoint)> functionToIntegrate = [phi, gradPhi, normal](RefPoint p) {
				return Utils::InnerProduct<Dim>(gradPhi(p), normal) * phi(p);
			};

			int polynomialDegree = reconstructPhi->GetDegree() - 1 + cellPhi->GetDegree();
			double integralFace = Kappa * face->ComputeIntegral(functionToIntegrate, polynomialDegree);

			sumFaces += integralFace;
		}

		return integralGradGrad - sumFaces;
	}

	double IntegrationByPartsRHS_face(Poisson_HHO_Face<Dim>* face, BasisFunction<Dim>* reconstructPhi, BasisFunction<Dim-1>* facePhi)
	{
		if (reconstructPhi->GetDegree() == 0)
			return 0;

		auto gradPhi = this->GradPhiOnFace(face, reconstructPhi);
		auto normal = this->OuterNormalVector(face);

		std::function<double(RefPoint)> functionToIntegrate = [facePhi, gradPhi, normal](RefPoint p) {
			return Utils::InnerProduct<Dim>(gradPhi(p), normal) * facePhi->Eval(p);
		};

		int polynomialDegree = reconstructPhi->GetDegree() - 1 + facePhi->GetDegree();
		return Kappa * face->ComputeIntegral(functionToIntegrate, polynomialDegree);
	}

	//---------------------------------------------//
	//--------------- Stabilization ---------------//
	//---------------------------------------------//

	void AssembleStabilizationMatrix()
	{
		this->Astab = Eigen::MatrixXd::Zero(CellBasis->Size() + this->Faces.size() * FaceBasis->Size(), CellBasis->Size() + this->Faces.size() * FaceBasis->Size());

		Eigen::MatrixXd ProjT = _projFromReconstruct;
		Eigen::MatrixXd Dt = ProjT * this->P;
		for (int i = 0; i < Dt.rows(); i++)
			Dt(i, i) -= 1;

		for (auto f : this->Faces)
		{
			Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
			Eigen::MatrixXd Mf = face->FaceMassMatrix();
			Eigen::MatrixXd ProjF = face->GetProjFromReconstruct(this);
			Eigen::MatrixXd ProjFT = face->GetProjFromCell(this);

			Eigen::MatrixXd Df = ProjF * this->P;
			for (int i = 0; i < Df.rows(); i++)
				Df(i, FirstDOFNumber(face) + i) -= 1;

			Eigen::MatrixXd DiffTF = Df - ProjFT * Dt;
			double h = face->GetDiameter();
			this->Astab += DiffTF.transpose() * Mf * DiffTF * Kappa / h;
		}
	}
	int DOFNumber(BasisFunction<Dim> * cellPhi)
	{
		return cellPhi->LocalNumber;
	}
	int DOFNumber(Face<Dim> * face, BasisFunction<Dim - 1> * facePhi)
	{
		return FirstDOFNumber(face) + facePhi->LocalNumber;
	}

public:
	int FirstDOFNumber(Face<Dim> * face)
	{
		return this->CellBasis->Size() + this->LocalNumberOf(face) * this->FaceBasis->Size();
	}
};