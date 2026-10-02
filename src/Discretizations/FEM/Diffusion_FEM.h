#pragma once
#include "../../Mesh/Mesh.h"
#include "../../Utils/Utils.h"
#include "../../TestCases/Diffusion/DiffusionTestCase.h"
#include "../../Utils/Parallelism.h"
#include "../../Utils/ExportModule.h"
#ifdef ENABLE_3D
#include "../../Geometry/3D/Tetrahedron.h"
#endif // ENABLE_3D
#include "../../Geometry/CartesianShape.h"
#include "../../Geometry/2D/Triangle.h"
using namespace std;

template <int Dim>
class Diffusion_FEM
{
private:
	FunctionalBasis<Dim>* _basis;
public:
	Mesh<Dim>* _mesh;
	DiffusionTestCase<Dim>* TestCase;
	
	Vector x_d; // Solution on the Dirichlet faces

	SparseMatrix A; // = A_ii
	SparseMatrix A_id;
	Vector b;

public:
	Diffusion_FEM() {}

	Diffusion_FEM(Mesh<Dim>* mesh, DiffusionTestCase<Dim>* testCase, FunctionalBasis<Dim>* basis)
	{	
		this->_mesh = mesh;
		this->TestCase = testCase;
		this->_basis = basis;

		// Re-numbering of the vertices: interior first, then Neumann, and Dirichlet at the end (because they will be eliminated from the system)
		BigNumber number = 0;
		for (auto v : this->_mesh->InteriorVertices)
			v->Number = number++;
		for (auto v : this->_mesh->NeumannVertices)
			v->Number = number++;
		for (auto v : this->_mesh->DirichletVertices)
			v->Number = number++;
	}

	void PrintDiscretization()
	{
		cout << "Mesh: " << _mesh->Description() << endl;
		cout << "    Elements  : " << _mesh->Elements.size() << endl;
		cout << "    Vertices  : " << _mesh->Vertices.size() << " (" << _mesh->InteriorVertices.size() << " interior + " << _mesh->BoundaryVertices.size() << " boundary)" << endl;
		cout << "    h         : " << scientific << this->_mesh->H() << defaultfloat << endl;
		cout << "    Regularity: " << this->_mesh->Regularity() << defaultfloat << endl;
		cout << "Discretization: continuous FEM (p = 1)" << endl;
		cout << "Unknowns : " << (_mesh->InteriorVertices.size() + _mesh->NeumannVertices.size()) << endl;
	}

	//--------------------------------------------//
	// Performs the assembly of the linear system //
	//--------------------------------------------//

	void Assemble(const ActionsArguments& actions)
	{
		Assemble(actions, ExportModule(""));
	}

	void Assemble(const ActionsArguments& actions, const ExportModule& out)
	{
		if (actions.LogAssembly)
			this->PrintDiscretization();

		if (actions.LogAssembly)
			cout << endl << "Assembly..." << endl;

		// Compute integrals on reference elements //
		if (actions.InitReferenceShapes)
			InitReferenceShapes();

		// Parallel loop on the elements //
		ThreadLocalCoeffs coeffs;
		auto basisFunctions = _basis->LocalFunctions();
		#pragma omp parallel for
		for (Element<Dim>* e : _mesh->Elements)
		{
			vector<Vertex*> vertices = e->Vertices();

			for (int i = 0; i< vertices.size(); i++)
			{
				Vertex* vi = vertices[i];
				auto phi_i = basisFunctions[i];
				for (int j = i; j < vertices.size(); j++)
				{
					Vertex* vj = vertices[j];
					auto phi_j = basisFunctions[j];

					double coeff = e->IntegralKGradGrad(e->DiffTensor(), phi_i, phi_j);
					coeffs.Local().Add(vi->Number, vj->Number, coeff);
					if (j != i)
						coeffs.Local().Add(vj->Number, vi->Number, coeff);
				}
			}
		}

		SparseMatrix globalA(_mesh->Vertices.size(), _mesh->Vertices.size());
		coeffs.Fill(globalA);

		// A = A_ii
		this->A    = globalA.topLeftCorner(_mesh->InteriorVertices.size(), _mesh->InteriorVertices.size());
		this->A_id = globalA.topRightCorner(_mesh->InteriorVertices.size(), _mesh->DirichletVertices.size());


		// Right-hand side //
		if (actions.AssembleRightHandSide)
		{
			Vector b_i = AssembleSourceTerm(TestCase->SourceFunction);
			this->x_d  = std::move(AssembleDirichletTerm(TestCase->BC.DirichletFunction));

			this->b = std::move(RHS(b_i, x_d)); //b_i - A_id * x_d;
		}
	}

	void InitReferenceShapes()
	{
		InitReferenceShapes(&this->TestCase->DiffField);
	}

	// Compute some useful integrals on reference element and store them
	void InitReferenceShapes(DiffusionField<Dim>* diffField)
	{
		// - Cartesian element
		CartesianShape<Dim, Dim>::InitReferenceShape()->ComputeAndStoreMassMatrix(_basis);
		if (Dim == 2)
		{
			// - Triangular element
			Triangle::InitReferenceShape()->ComputeAndStoreMassMatrix((FunctionalBasis<2>*)_basis);
		}
		else if (Dim == 3)
		{
#ifdef ENABLE_3D
			// - Tetrahedral element
			Tetrahedron::InitReferenceShape()->ComputeAndStoreMassMatrix((FunctionalBasis<3>*)_basis);
#endif // ENABLE_3D
		}
	}

	SparseMatrix MassMatrix()
	{
		ThreadLocalCoeffs coeffs;
		#pragma omp parallel for
		for (Element<Dim>* e : _mesh->Elements)
		{
			DenseMatrix mass = e->MassMatrix(_basis);
			vector<Vertex*> vertices = e->Vertices();

			for (int i = 0; i < vertices.size(); i++)
			{
				Vertex* vi = vertices[i];
				//auto phi_i = _basis->LocalFunctions[i];
				for (int j = i; j < vertices.size(); j++)
				{
					Vertex* vj = vertices[j];
					//auto phi_j = _basis->LocalFunctions[j];

					double coeff = mass(i, j);
					coeffs.Local().Add(vi->Number, vj->Number, coeff);
					if (j != i)
						coeffs.Local().Add(vj->Number, vi->Number, coeff);
				}
			}
		}

		SparseMatrix mass(_mesh->Vertices.size(), _mesh->Vertices.size());
		coeffs.Fill(mass);
		return mass;
	}

	//------------------------------//
	//   Assemble right-hand side   //
	//------------------------------//

	Vector AssembleSourceTerm(DomFunction sourceFunction)
	{
		Vector b_i = Vector::Zero(_mesh->InteriorVertices.size());
		auto localFunctions = _basis->LocalFunctions();
		#pragma omp parallel for
		for (Element<Dim>* e : this->_mesh->Elements)
		{
			vector<Vertex*> vertices = e->Vertices();
			for (int i = 0; i < vertices.size(); i++)
			{
				Vertex* v = vertices[i];
				if (v->Number < _mesh->InteriorVertices.size())
				{
					double value = e->SourceTerm(localFunctions[i], sourceFunction);
					#pragma omp atomic // the vertices are shared by several elements
					b_i[v->Number] += value;
				}
			}
		}

		return b_i;
	}

	Vector AssembleSourceTerm(const Vector& sourceNodeValues, const SparseMatrix& massMatrix)
	{
		return massMatrix.topLeftCorner(_mesh->InteriorVertices.size(), _mesh->InteriorVertices.size()) * sourceNodeValues.head(_mesh->InteriorVertices.size());
	}

	// Solution on the Dirichlet nodes
	Vector AssembleDirichletTerm(DomFunction dirichletFunction)
	{
		Vector x_d = Vector(_mesh->DirichletVertices.size());
		#pragma omp parallel for
		for (Vertex* v : this->_mesh->DirichletVertices)
			x_d[v->Number - _mesh->InteriorVertices.size()] = dirichletFunction(*v);
		return x_d;
	}

	Vector RHS(const Vector& b_i, const Vector& x_d)
	{
		return b_i - A_id * x_d;
	}

	//------------------------------//
	//        Final solution        //
	//------------------------------//

	Vector BuildCompleteSolution(const Vector& x_i, const Vector& x_d)
	{
		Vector x(x_i.rows() + x_d.rows());
		x << x_i, x_d;
		return x;
	}

	Vector AddDirichletValues(const Vector& x_i)
	{
		return BuildCompleteSolution(x_i, this->x_d);
	}

	//------------------------------//
	//           L2-error           //
	//------------------------------//

	double L2Error(DomFunction exactSolution, const Vector& solution)
	{
		assert(solution.rows() == _mesh->Vertices.size());

		double absoluteError = 0;
		double normExactSolution = 0;

		auto localFunctions = _basis->LocalFunctions();
		#pragma omp parallel for reduction(+:absoluteError, normExactSolution)
		for (Element<Dim>* e : this->_mesh->Elements)
		{
			RefFunction approximate = [this, &solution, &localFunctions, e](const RefPoint& p) {
				double total = 0;
				vector<Vertex*> vertices = e->Vertices();
				for (int i = 0; i < vertices.size(); i++)
				{
					Vertex* v = vertices[i];
					total += solution(v->Number) * localFunctions[i]->Eval(p);
				}
				return total;
			};

			absoluteError += e->L2ErrorPow2(approximate, exactSolution);
			normExactSolution += e->Integral([exactSolution](const DomPoint& p) { return pow(exactSolution(p), 2); });
		}

		absoluteError = sqrt(absoluteError);
		normExactSolution = sqrt(normExactSolution);
		return normExactSolution != 0 ? absoluteError / normExactSolution : absoluteError;
	}

	void AssertSchemeConvergence(double l2Error)
	{
		if (Dim == 1)
			return;

		double h = this->_mesh->H();
		double constant = 100;
		double upperBound;
		if ((TestCase->Code().compare("sine") == 0
			|| TestCase->Code().compare("poly") == 0
			|| TestCase->Code().compare("zero") == 0) && TestCase->DiffField.IsHomogeneous)
		{
			upperBound = constant * pow(h, 2);
		}
		else
			return;

		cout << "UpperBound = " << upperBound << endl;
		if (l2Error > upperBound)
			cout << Utils::BeginRed << "The L2 error is bigger than the theoretical upper bound." << Utils::EndColor << endl;
	}

/*public:

	double IntegralOverDomain(DomFunction func)
	{
		double total = 0;
		#pragma omp parallel for reduction(+:total)
		for (Element<Dim>* e : _mesh->Elements)
			total += e->Integral(func);
		return total;
	}

	double IntegralOverBoundary(DomFunction func)
	{
		double total = 0;
		#pragma omp parallel for reduction(+:total)
		for (Face<Dim>* f : _mesh->BoundaryFaces)
			total += f->Integral(func);
		return total;
	}*/
};

