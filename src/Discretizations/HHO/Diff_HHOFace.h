#pragma once
#include "../../FunctionalBasis/Orthogonal/OrthogonalBasisOnCstJacShape.h"
#include "HHOParameters.h"

template <int Dim>
class Diff_HHOElement;

template <int Dim>
class Diff_HHOFace
{
private:
	DenseMatrix _massMatrix;
	Eigen::LLT<DenseMatrix> _massMatrixSolver;
public:
	Face<Dim>* MeshFace = nullptr;
	HHOParameters<Dim>* HHO;
	FunctionalBasis<Dim - 1>* Basis = nullptr;

	Diff_HHOFace() {}

	//--------------//
	// Face wrapper //
	//--------------//

public:
	inline BigNumber Number() const
	{
		return this->MeshFace->Number;
	}
	inline double Diameter() const
	{
		return this->MeshFace->Diameter();
	}
	inline bool HasDirichletBC()
	{
		return this->MeshFace->HasDirichletBC();
	}
	inline bool HasNeumannBC()
	{
		return this->MeshFace->HasNeumannBC();
	}
private:
	inline DomPoint ConvertToDomain(const RefPoint& refPoint) const
	{
		return this->MeshFace->ConvertToDomain(refPoint);
	}
	inline RefPoint ConvertToReference(const DomPoint& domainPoint) const
	{
		return this->MeshFace->ConvertToReference(domainPoint);
	}

	//-----------//
	//    HHO    //
	//-----------//

public:
	void InitHHO(HHOParameters<Dim>* hho)
	{
		if (this->Basis)
			return;

		this->HHO = hho;

		//this->ComputeAndSaveQuadraturePoints(basis->GetDegree());
		//this->ComputeAndSaveQuadraturePoints();

		if (hho->OrthogonalizeFaceBases())
		{
			this->Basis = new OrthogonalBasis<Dim - 1>(HHO->FaceBasis, this->MeshFace->Shape(), hho->NFaceOrthogonalizations(), hho->OrthonormalizeFaceBases());
			//this->_massMatrix = this->MeshFace->Shape()->MassMatrix(this->Basis);
			//cout << "mass matrix: " << endl << _massMatrix << endl;
		}
		else
		{
			this->Basis = hho->FaceBasis;
			this->_massMatrix = this->MeshFace->Shape()->MassMatrix(this->Basis);
			this->_massMatrixSolver = this->_massMatrix.llt();
		}
	}

	DenseMatrix MassMatrix()
	{
		if (HHO->OrthonormalizeFaceBases())
			return DenseMatrix::Identity(this->Basis->Size(), this->Basis->Size());
		else if (HHO->OrthogonalizeFaceBases())
		{
			Vector d(this->Basis->Size());
			for (BasisFunction<Dim-1>* phi : this->Basis->LocalFunctions())
				d[phi->LocalNumber] = dynamic_cast<OrthogonalBasisFunction<Dim-1>*>(phi)->NormSquare();
			return d.asDiagonal();
		}
		else
			return _massMatrix;
	}
	Vector ApplyMassMatrix(const Vector& v)
	{
		if (HHO->OrthonormalizeFaceBases())
			return v;
		else if (HHO->OrthogonalizeFaceBases())
		{
			Vector d(this->Basis->Size());
			for (BasisFunction<Dim-1>* phi : this->Basis->LocalFunctions())
				d[phi->LocalNumber] = dynamic_cast<OrthogonalBasisFunction<Dim-1>*>(phi)->NormSquare();
			return d.asDiagonal() * v;
		}
		else
		{
			assert(_massMatrix.rows() > 0 && "The face mass matrix has not been initialized or was deleted.");
			return _massMatrix * v;
		}
	}
	DenseMatrix SolveMassMatrix(const DenseMatrix& M) // also works for vectors
	{
		if (HHO->OrthonormalizeFaceBases())
			return M;
		else if (HHO->OrthogonalizeFaceBases())
		{
			Vector d(this->Basis->Size());
			for (BasisFunction<Dim-1>* phi : this->Basis->LocalFunctions())
				d[phi->LocalNumber] = dynamic_cast<OrthogonalBasisFunction<Dim-1>*>(phi)->NormSquare();
			return d.asDiagonal().inverse() * M;
		}
		else
			return _massMatrixSolver.solve(M);
	}
	// Included in the preceding function!
	/*Vector SolveMassMatrix(const Vector& v)
	{
		if (HHO->OrthonormalizeFaceBases())
			return v;
		else if (HHO->OrthogonalizeFaceBases())
		{
			Vector d(this->Basis->Size());
			for (BasisFunction<Dim-1>* phi : this->Basis->LocalFunctions)
				d[phi->LocalNumber] = dynamic_cast<OrthogonalBasisFunction<Dim-1>*>(phi)->NormSquare();
			return d.asDiagonal().inverse() * v;
		}
		else
			return _massMatrixSolver.solve(v);
	}*/

	double ComputeMassTerm(BasisFunction<Dim - 1>* facePhi, Element<Dim>* element, BasisFunction<Dim>* cellPhi)
	{
		RefFunction functionToIntegrate = [this, element, facePhi, cellPhi](const RefPoint& p) {
			return facePhi->Eval(p) * element->EvalTrace(this->MeshFace, cellPhi, p);
		};

		int polynomialDegree = facePhi->GetDegree() + cellPhi->GetDegree();
		return this->MeshFace->Integral(functionToIntegrate, polynomialDegree);
	}

private:
	double NormalDerivativeTerm(BasisFunction<Dim - 1>* facePhi, Element<Dim>* element, BasisFunction<Dim>* cellPhi, const DimVector<Dim>& n)
	{
		if (cellPhi->GetDegree() == 0)
			return 0;

		RefFunction functionToIntegrate = [this, element, facePhi, cellPhi, &n](const RefPoint& p) {
			return facePhi->Eval(p) * element->EvalGradOnFace(this->MeshFace, cellPhi, p).dot(n);
		};

		int polynomialDegree = facePhi->GetDegree() + cellPhi->GetDegree() - 1;
		return this->MeshFace->Integral(functionToIntegrate, polynomialDegree);
	}

public:
	DenseMatrix MassMatrix(FunctionalBasis<Dim - 1>* other)
	{
		DenseMatrix M(this->Basis->Size(), other->Size());
		for (BasisFunction<Dim - 1>* phi1 : this->Basis->LocalFunctions())
		{
			for (BasisFunction<Dim - 1>* phi2 : other->LocalFunctions())
			{
				double term = this->MeshFace->Shape()->ComputeMassTerm(phi1, phi2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
			}
		}
		return M;
	}

	DenseMatrix Trace(Element<Dim>* element, FunctionalBasis<Dim>* cellBasis)
	{
		DenseMatrix massFaceCell = this->MassMatrix(this->Basis, element, cellBasis);
		return SolveMassMatrix(massFaceCell);
	}

	DenseMatrix NormalDerivative(Element<Dim>* element, FunctionalBasis<Dim>* cellBasis)
	{
		DenseMatrix M(this->Basis->Size(), cellBasis->Size());

		DimVector<Dim> n = element->OuterNormalVector(this->MeshFace);
		for (BasisFunction<Dim-1>* facePhi : this->Basis->LocalFunctions())
		{
			for (BasisFunction<Dim>* cellPhi : cellBasis->LocalFunctions())
			{
				double term = this->NormalDerivativeTerm(facePhi, element, cellPhi, n);
				M(facePhi->LocalNumber, cellPhi->LocalNumber) = term;
			}
		}
		return SolveMassMatrix(M);
	}

	double Integral(const Vector& coeffs)
	{
		double integral = 0;
		for (BasisFunction<Dim-1>* phi : this->Basis->LocalFunctions())
			integral += coeffs[phi->LocalNumber] * this->MeshFace->Integral(phi); // TODO: this can be computed only once on the reference element
		return integral;
	}

private:
	double InnerProduct(BasisFunction<Dim - 1>* phi, DomFunction f)
	{
		RefFunction functionToIntegrate = [this, f, phi](const RefPoint& refElementPoint) {
			DomPoint domainPoint = this->ConvertToDomain(refElementPoint);
			return f(domainPoint) * phi->Eval(refElementPoint);
		};

		return this->MeshFace->Integral(functionToIntegrate);
	}

public:
	Vector InnerProductWithBasis(DomFunction f)
	{
		Vector innerProducts(this->Basis->Size());
		for (BasisFunction<Dim - 1>* phi : this->Basis->LocalFunctions())
			innerProducts(phi->LocalNumber) = InnerProduct(phi, f);
		return innerProducts;
	}
	Vector ProjectOnBasis(DomFunction f)
	{
		return SolveMassMatrix(InnerProductWithBasis(f));
	}

	DenseMatrix ProjectOnBasis(FunctionalBasis<Dim - 1>* other)
	{
		return SolveMassMatrix(MassMatrix(other));
	}

	double InnerProd(const Vector& coeffs1, const Vector& coeffs2)
	{
		assert(coeffs1.rows() == HHO->nFaceUnknowns);
		return coeffs1.dot(ApplyMassMatrix(coeffs2));
	}

	void DeleteUselessMatricesAfterAssembly()
	{
		if (!HHO->OrthogonalizeFaceBases())
			Utils::Empty(_massMatrix);
	}

	void DeleteUselessMatricesAfterMultigridSetup()
	{
		//if (!HHO->OrthogonalizeFaceBases())   // I comment it quick and dirty because I need this massMatrixSolver for the bi-harmonic problem
			//_massMatrixSolver = Eigen::LLT<DenseMatrix>();
		this->MeshFace->EmptySavedDomPoints();
	}

private:
	DenseMatrix MassMatrix(FunctionalBasis<Dim - 1>* basis, Element<Dim>* element, FunctionalBasis<Dim>* cellBasis)
	{
		DenseMatrix M(basis->Size(), cellBasis->Size());
		for (BasisFunction<Dim - 1>* phi1 : basis->LocalFunctions())
		{
			for (BasisFunction<Dim>* phi2 : cellBasis->LocalFunctions())
			{
				double term = this->ComputeMassTerm(phi1, element, phi2);
				M(phi1->LocalNumber, phi2->LocalNumber) = term;
			}
		}
		return M;
	}

public:
	~Diff_HHOFace()
	{
		if (HHO->OrthogonalizeFaceBases())
			delete Basis;
	}
};