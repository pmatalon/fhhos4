#pragma once
#include <Eigen/Sparse>
#include "../Mesh/Element.h"
#include "../Utils/SourceFunction.h"
#include "HHOParameters.h"
#include "Poisson_HHO_Face.h"

template <int Dim>
class Poisson_HHO_Element : virtual public Element<Dim>
{
private:
	Eigen::MatrixXd _cellMassMatrix;
	Eigen::MatrixXd _projFromReconstruct;
public:
	HHOParameters<Dim>* HHO;

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
	
	void InitHHO(HHOParameters<Dim>* hho)
	{
		this->HHO = hho;

		this->_cellMassMatrix = this->CellMassMatrix(hho->CellBasis);
		Eigen::MatrixXd Nt = this->CellReconstructMassMatrix(hho->CellBasis, hho->ReconstructionBasis);
		this->_projFromReconstruct = _cellMassMatrix.inverse() * Nt;

		this->AssembleReconstructionAndConsistencyMatrices();
		this->AssembleStabilizationMatrix();

		int nTotalFaceUnknowns = this->Faces.size() * hho->nFaceUnknowns;

		this->A = Acons + Astab;
		auto Att = A.topLeftCorner(hho->nCellUnknowns, hho->nCellUnknowns);
		//auto Aff = A.bottomRightCorner(nTotalFaceUnknowns, nTotalFaceUnknowns);
		//auto Atf = A.topRightCorner(nCellUnknowns, nTotalFaceUnknowns);
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
		int nTotalFaceUnknowns = this->Faces.size() * HHO->nFaceUnknowns;

		auto Atf = this->A.topRightCorner(HHO->nCellUnknowns, nTotalFaceUnknowns);
		Eigen::MatrixXd solveCellUnknowns = -this->invAtt * Atf;

		Eigen::MatrixXd createHybridVectorFromFacesMatrix(HHO->nCellUnknowns + nTotalFaceUnknowns, nTotalFaceUnknowns);
		createHybridVectorFromFacesMatrix.topRows(HHO->nCellUnknowns) = solveCellUnknowns;
		createHybridVectorFromFacesMatrix.bottomRows(nTotalFaceUnknowns) = Eigen::MatrixXd::Identity(nTotalFaceUnknowns, nTotalFaceUnknowns);
		return this->P * createHybridVectorFromFacesMatrix;
	}

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
	virtual double IntegralKGradGradReconstruct(Tensor<Dim>* K, BasisFunction<Dim>* reconstructPhi1, BasisFunction<Dim>* reconstructPhi2) = 0;
	virtual double ComputeIntegralKGradGrad(Tensor<Dim>* K, BasisFunction<Dim>* phi1, BasisFunction<Dim>* phi2) const = 0;
	virtual Eigen::MatrixXd CellMassMatrix(FunctionalBasis<Dim>* basis) = 0;
	virtual Eigen::MatrixXd CellReconstructMassMatrix(FunctionalBasis<Dim>* cellBasis, FunctionalBasis<Dim>* reconstructBasis) = 0;


	virtual ~Poisson_HHO_Element() {}

private:

	//------------------------------------------------------------------------------------//
	//--------------- Reconstruction operator and consistency contribution ---------------//
	//------------------------------------------------------------------------------------//

	void AssembleReconstructionAndConsistencyMatrices()
	{
		this->P = Eigen::MatrixXd(HHO->nReconstructUnknowns, HHO->nCellUnknowns + this->Faces.size() * HHO->nFaceUnknowns);

		Eigen::MatrixXd reconstructionMatrixToInvert = this->AssembleReconstructionMatrixToInvert();
		Eigen::MatrixXd rhsMatrix = this->AssembleRHSMatrix();

		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> solver = reconstructionMatrixToInvert.colPivHouseholderQr();
		for (int j = 0; j < rhsMatrix.cols(); j++)
		{
			Eigen::VectorXd col = solver.solve(rhsMatrix.col(j));

			// We don't keep the last element of the column, which is the Lagrange multiplier
			auto colWithoutLagrangeMultiplier = col.head(HHO->nReconstructUnknowns);
			this->P.col(j) = colWithoutLagrangeMultiplier;
		}

		Eigen::MatrixXd laplacianMatrix = reconstructionMatrixToInvert.topLeftCorner(HHO->nReconstructUnknowns, HHO->nReconstructUnknowns); // Grad-Grad part (block S)

		this->Acons = this->P.transpose() * laplacianMatrix * this->P;
	}

	Eigen::MatrixXd AssembleReconstructionMatrixToInvert()
	{
		Eigen::MatrixXd matrixToInvert(HHO->nReconstructUnknowns + 1, HHO->nReconstructUnknowns + 1);
		this->AssembleGradientReconstructionMatrix(matrixToInvert); // Block S
		this->AssembleMeanValueCondition(matrixToInvert); // Blocks L and L_transpose
		matrixToInvert.bottomRightCorner<1, 1>() << 0;
		return matrixToInvert;
	}

	void AssembleGradientReconstructionMatrix(Eigen::MatrixXd & reconstructionMatrixToInvert)
	{
		for (BasisFunction<Dim>* phi1 : HHO->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim>* phi2 : HHO->ReconstructionBasis->LocalFunctions)
				reconstructionMatrixToInvert(phi1->LocalNumber, phi2->LocalNumber) = this->IntegralKGradGradReconstruct(this->DiffTensor, phi1, phi2);
		}
	}

	void AssembleMeanValueCondition(Eigen::MatrixXd & reconstructionMatrixToInvert)
	{
		int last = HHO->ReconstructionBasis->NumberOfLocalFunctionsInElement(NULL);
		for (BasisFunction<Dim>* phi : HHO->ReconstructionBasis->LocalFunctions)
		{
			double meanValue = this->Integral(phi);
			reconstructionMatrixToInvert(last, phi->LocalNumber) = meanValue;
			reconstructionMatrixToInvert(phi->LocalNumber, last) = meanValue;
		}
	}

	void AssembleMeanValueConditionRHS(Eigen::MatrixXd & rhsMatrix)
	{
		int last = HHO->ReconstructionBasis->NumberOfLocalFunctionsInElement(NULL);
		for (BasisFunction<Dim>* phi : HHO->CellBasis->LocalFunctions)
			rhsMatrix(last, phi->LocalNumber) = this->Integral(phi);
	}

	Eigen::MatrixXd AssembleRHSMatrix()
	{
		auto nTotalFaceUnknowns = this->Faces.size() * HHO->nFaceUnknowns;
		auto nColumns = HHO->nCellUnknowns + nTotalFaceUnknowns;
		Eigen::MatrixXd rhsMatrix(HHO->nReconstructUnknowns + 1, nColumns);

		// Top-left corner (Block Bt)
		this->AssembleIntegrationByPartsRHS_cell(rhsMatrix);

		// Top-right corner (Block B_frontier)
		for (auto face : this->Faces)
			this->AssembleIntegrationByPartsRHS_face(rhsMatrix, face);

		// Bottom-left corner (Block Lt)
		this->AssembleMeanValueConditionRHS(rhsMatrix);

		// Bottom-right corner (0)
		rhsMatrix.bottomRightCorner(1, nTotalFaceUnknowns) << Eigen::ArrayXXd::Zero(1, nTotalFaceUnknowns);

		return rhsMatrix;
	}

	void AssembleIntegrationByPartsRHS_cell(Eigen::MatrixXd & rhsMatrix)
	{
		for (BasisFunction<Dim>* reconstructPhi : HHO->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim>* cellPhi : HHO->CellBasis->LocalFunctions)
				rhsMatrix(reconstructPhi->LocalNumber, DOFNumber(cellPhi)) = this->IntegrationByPartsRHS_cell(reconstructPhi, cellPhi);
		}
	}

	void AssembleIntegrationByPartsRHS_face(Eigen::MatrixXd & rhsMatrix, Face<Dim> * f)
	{
		Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
		for (BasisFunction<Dim>* reconstructPhi : HHO->ReconstructionBasis->LocalFunctions)
		{
			for (BasisFunction<Dim - 1> * facePhi : HHO->FaceBasis->LocalFunctions)
				rhsMatrix(reconstructPhi->LocalNumber, DOFNumber(face, facePhi)) = this->IntegrationByPartsRHS_face(face, reconstructPhi, facePhi);
		}
	}

	double IntegrationByPartsRHS_cell(BasisFunction<Dim>* reconstructPhi, BasisFunction<Dim>* cellPhi)
	{
		if (reconstructPhi->GetDegree() == 0)
			return 0;

		double integralGradGrad = this->ComputeIntegralKGradGrad(this->DiffTensor, reconstructPhi, cellPhi);

		double sumFaces = 0;
		for (auto f : this->Faces)
		{
			Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);

			auto phi = this->EvalPhiOnFace(face, cellPhi);
			auto gradPhi = this->GradPhiOnFace(face, reconstructPhi);
			auto normal = this->OuterNormalVector(face);

			std::function<double(RefPoint)> functionToIntegrate = [this, phi, gradPhi, normal](RefPoint p) {
				return (this->DiffTensor * gradPhi(p)).dot(normal) * phi(p);
			};

			int polynomialDegree = reconstructPhi->GetDegree() - 1 + cellPhi->GetDegree();
			double integralFace = face->ComputeIntegral(functionToIntegrate, polynomialDegree);

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

		std::function<double(RefPoint)> functionToIntegrate = [this, facePhi, gradPhi, normal](RefPoint p) {
			return (this->DiffTensor * gradPhi(p)).dot(normal) * facePhi->Eval(p);
		};

		int polynomialDegree = reconstructPhi->GetDegree() - 1 + facePhi->GetDegree();
		return face->ComputeIntegral(functionToIntegrate, polynomialDegree);
	}

	//---------------------------------------------//
	//--------------- Stabilization ---------------//
	//---------------------------------------------//

	void AssembleStabilizationMatrix()
	{
		int nHybridUnknowns = HHO->nCellUnknowns + this->Faces.size() * HHO->nFaceUnknowns;

		this->Astab = Eigen::MatrixXd::Zero(nHybridUnknowns, nHybridUnknowns);

		if (HHO->Stabilization.compare("hho") == 0)
		{
			Eigen::MatrixXd ProjT = _projFromReconstruct;
			Eigen::MatrixXd Dt = ProjT * this->P;
			for (int i = 0; i < Dt.rows(); i++)
				Dt(i, i) -= 1;

			for (auto f : this->Faces)
			{
				Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
				auto normal = this->OuterNormalVector(face);
				Eigen::MatrixXd Mf = face->FaceMassMatrix();
				Eigen::MatrixXd ProjF = face->GetProjFromReconstruct(this);
				Eigen::MatrixXd ProjFT = face->GetProjFromCell(this);

				Eigen::MatrixXd Df = ProjF * this->P;
				for (int i = 0; i < Df.rows(); i++)
					Df(i, FirstDOFNumber(face) + i) -= 1;

				Eigen::MatrixXd DiffTF = Df - ProjFT * Dt;
				double h = face->GetDiameter();
				this->Astab += DiffTF.transpose() * Mf * DiffTF * (this->DiffTensor * normal).dot(normal) / h;
			}
		}
		else if (HHO->Stabilization.compare("hdg") == 0)
		{
			for (auto f : this->Faces)
			{
				Poisson_HHO_Face<Dim>* face = dynamic_cast<Poisson_HHO_Face<Dim>*>(f);
				auto normal = this->OuterNormalVector(face);
				Eigen::MatrixXd Mf = face->FaceMassMatrix();
				Eigen::MatrixXd ProjFT = face->GetProjFromCell(this);

				Eigen::MatrixXd Fpart = Eigen::MatrixXd::Zero(HHO->nFaceUnknowns, nHybridUnknowns);
				Fpart.middleCols(FirstDOFNumber(face), HHO->nFaceUnknowns) = Eigen::MatrixXd::Identity(HHO->nFaceUnknowns, HHO->nFaceUnknowns);

				Eigen::MatrixXd Tpart = Eigen::MatrixXd::Zero(HHO->nCellUnknowns, nHybridUnknowns);
				Tpart.middleCols(0, HHO->nCellUnknowns) = Eigen::MatrixXd::Identity(HHO->nCellUnknowns, HHO->nCellUnknowns);

				Eigen::MatrixXd DiffTF = Fpart - ProjFT * Tpart;
				double h = face->GetDiameter();
				this->Astab += DiffTF.transpose() * Mf * DiffTF * (this->DiffTensor * normal).dot(normal) / h;
			}
		}
		else
			assert(false);
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
		return HHO->nCellUnknowns + this->LocalNumberOf(face) * HHO->nFaceUnknowns;
	}
};