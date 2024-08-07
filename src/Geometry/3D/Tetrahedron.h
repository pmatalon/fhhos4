#pragma once
#include "../../Mesh/Vertex.h"
#include "ReferenceTetrahedron.h"
#include "../PhysicalShapeWithConstantJacobian.h"
using namespace std;

class Tetrahedron : public PhysicalShapeWithConstantJacobian<3>
{
private:
	DomPoint v1;
	DomPoint v2;
	DomPoint v3;
	DomPoint v4;

	double _diameter;
	double _measure;
	DomPoint _center;
	double _inRadius;

	DimMatrix<3> _inverseMapping;
	DimMatrix<3> _inverseJacobianTranspose;
	double _detJacobian;

	vector<Tetrahedron> _refinement;

public:
	static ReferenceTetrahedron RefTetra;

	Tetrahedron() {}

	Tetrahedron(const DomPoint& p1, const DomPoint& p2, const DomPoint& p3, const DomPoint& p4)
		: v1(p1), v2(p2), v3(p3), v4(p4), _refinement(0)
	{
		assert(p1 != p2 && p1 != p3 && p1 != p4 && p2 != p3 && p2 != p4 && p3 != p4);
		Init();
	}

	Tetrahedron(const Tetrahedron& shape) = default;

	inline DomPoint V1() const { return v1; }
	inline DomPoint V2() const { return v2; }
	inline DomPoint V3() const { return v3; }
	inline DomPoint V4() const { return v4; }

	void Init()
	{
		DimVector<3> v12 = Vect<3>(v1, v2);
		DimVector<3> v13 = Vect<3>(v1, v3);
		DimVector<3> v14 = Vect<3>(v1, v4);
		DimVector<3> v23 = Vect<3>(v2, v3);
		DimVector<3> v24 = Vect<3>(v2, v4);
		DimVector<3> v34 = Vect<3>(v3, v4);
		double lengthEdge12 = v12.norm();
		double lengthEdge13 = v13.norm();
		double lengthEdge14 = v14.norm();
		double lengthEdge23 = v23.norm();
		double lengthEdge24 = v24.norm();
		double lengthEdge34 = v34.norm();
		_diameter = max({ lengthEdge12, lengthEdge13, lengthEdge14, lengthEdge23, lengthEdge24, lengthEdge34 });

		DimMatrix<3> m;
		m.col(0) = v12;
		m.col(1) = v13;
		m.col(2) = v14;
		_measure = abs(m.determinant()) / 6;

		_center = DomPoint((v1.X + v2.X + v3.X + v4.X) / 4, (v1.Y + v2.Y + v3.Y + v4.Y) / 4, (v1.Z + v2.Z + v3.Z + v4.Z) / 4);

		_inRadius = 3 * _measure / (0.5*v12.cross(v13).norm() + 0.5*v13.cross(v14).norm() + 0.5*v14.cross(v12).norm() + 0.5*v23.cross(v24).norm());

		_detJacobian = _measure / RefTetra.Measure();


		DimMatrix<3> mapping;
		mapping <<
			v2.X - v1.X, v3.X - v1.X, v4.X - v1.X,
			v2.Y - v1.Y, v3.Y - v1.Y, v4.Y - v1.Y,
			v2.Z - v1.Z, v3.Z - v1.Z, v4.Z - v1.Z;
		_inverseMapping = mapping.inverse();

		DimMatrix<3> inverseJacobian = _inverseMapping;
		_inverseJacobianTranspose = inverseJacobian.transpose();

		//assert(abs(_detJacobian - 1 / inverseJacobian.determinant()) < 1e-14);
	}

	PhysicalShape<3>* CreateCopy() const override
	{
		return new Tetrahedron(*this);
	}

	ReferenceShape<3>* RefShape() const override
	{
		return &RefTetra;
	}

	bool MapsToACartesianShape() const override
	{
		return false;
	}
	
	inline vector<DomPoint> Vertices() const override
	{
		return vector<DomPoint>{ v1, v2, v3, v4 };
	}

	bool IsDegenerated() const override
	{
		Utils::FatalError("The function Tetrahedron::IsDegenerated() is not implemented.");
		return true;
	}

	static ReferenceTetrahedron* InitReferenceShape()
	{
		return &RefTetra;
	}

	inline double Diameter() const override
	{
		return _diameter;
	}
	inline double Measure() const override
	{
		return _measure;
	}
	inline DomPoint Center() const override
	{
		return _center;
	}
	inline DomPoint InteriorPoint() const override
	{
		return _center;
	}
	inline bool IsConvex() const override
	{
		return true;
	}
	inline double InRadius() const override
	{
		return _inRadius;
	}
	inline bool Contains(const DomPoint& p) const override
	{
		return TetrahedronContains(v1, v2, v3, v4, p, _measure);
	}

	static bool TetrahedronContains(const DomPoint& A, const DomPoint& B, const DomPoint& C, const DomPoint& D, const DomPoint& P, double tetraVolume)
	{
		// Barycentric coordinates
		double alpha = Tetrahedron::Volume(P, B, C, D) / tetraVolume;
		double beta = Tetrahedron::Volume(A, P, C, D) / tetraVolume;
		double gamma = Tetrahedron::Volume(A, B, P, D) / tetraVolume;
		double delta = Tetrahedron::Volume(A, B, C, P) / tetraVolume;

		double tol = Utils::Eps;
		return alpha + tol > 0 && alpha < 1 + tol    // alpha >= 0 && alpha <= 1
			&& beta + tol > 0 && beta < 1 + tol    // beta  >= 0 && beta  <= 1
			&& gamma + tol > 0 && gamma < 1 + tol    // gamma >= 0 && gamma <= 1
			&& delta + tol > 0 && delta < 1 + tol    // delta >= 0 && delta <= 1
			&& (abs(alpha + beta + gamma + delta - 1) < tol); // alpha + beta + gamma + delta = 1
	}

	static double Volume(const DomPoint& A, const DomPoint& B, const DomPoint& C, const DomPoint& D)
	{
		DimMatrix<3> m;
		m.col(0) = Vect<3>(A, B);
		m.col(1) = Vect<3>(A, C);
		m.col(2) = Vect<3>(A, D);
		return abs(m.determinant()) / 6;
	}

	void RefineWithoutCoarseOverlap(const vector<PhysicalShape<2>*>& doNotCross) override
	{
		assert(false && "RefineWithoutCoarseOverlap() not implemented for the tetrahedron.");
		Utils::FatalError("RefineWithoutCoarseOverlap() not implemented for the tetrahedron.");
	}

	void RefineByBey(int nRefinements)
	{
		if (!_refinement.empty())
			return;

		vector<Tetrahedron> firstRefinement;
		firstRefinement.reserve(8);

		DomPoint m12 = Middle<3>(v1, v2);
		DomPoint m13 = Middle<3>(v1, v3);
		DomPoint m14 = Middle<3>(v1, v4);
		DomPoint m23 = Middle<3>(v2, v3);
		DomPoint m24 = Middle<3>(v2, v4);
		DomPoint m34 = Middle<3>(v3, v4);

		// Corners of the tetrahedron
		firstRefinement.emplace_back(v1, m12, m13, m14);
		firstRefinement.emplace_back(m12, v2, m23, m24);
		firstRefinement.emplace_back(m13, m23, v3, m34);
		firstRefinement.emplace_back(m14, m24, m34, v4);

		// Remaining octahedron
		firstRefinement.emplace_back(m12, m13, m14, m24);
		firstRefinement.emplace_back(m12, m13, m23, m24);
		firstRefinement.emplace_back(m13, m14, m24, m34);
		firstRefinement.emplace_back(m13, m23, m24, m34);

		if (nRefinements == 1)
			_refinement = firstRefinement;
		else
		{
			for (Tetrahedron& t1 : firstRefinement)
			{
				t1.RefineByBey(nRefinements - 1);
				for (Tetrahedron& t2 : t1._refinement)
					_refinement.push_back(t2);
			}
		}
	}

	vector<const PhysicalShape<3>*> RefinedShapes() const override
	{
		assert(_refinement.size() > 0);
		vector<const PhysicalShape<3>*> refinedShapes;
		for (const Tetrahedron& t : _refinement)
		{
			const PhysicalShape<3>* ps = &t;
			refinedShapes.push_back(ps);
		}
		return refinedShapes;
	}
	vector<PhysicalShape<3>*> RefinedShapes() override
	{
		assert(_refinement.size() > 0);
		vector<PhysicalShape<3>*> refinedShapes;
		for (Tetrahedron& t : _refinement)
		{
			PhysicalShape<3>* ps = &t;
			refinedShapes.push_back(ps);
		}
		return refinedShapes;
	}

	inline double DetJacobian() const override
	{
		return _detJacobian;
	}
	DimMatrix<3> InverseJacobianTranspose() const override
	{
		return _inverseJacobianTranspose;
	}

	// Mapping
	DomPoint ConvertToDomain(const RefPoint& refPoint) const override
	{
		double t = refPoint.X;
		double u = refPoint.Y;
		double v = refPoint.Z;

		DomPoint p;
		p.X = (v2.X - v1.X)*t + (v3.X - v1.X)*u + (v4.X - v1.X)*v + v1.X;
		p.Y = (v2.Y - v1.Y)*t + (v3.Y - v1.Y)*u + (v4.Y - v1.Y)*v + v1.Y;
		p.Z = (v2.Z - v1.Z)*t + (v3.Z - v1.Z)*u + (v4.Z - v1.Z)*v + v1.Z;
		return p;
	}

	// Inverse mapping
	RefPoint ConvertToReference(const DomPoint& domainPoint) const override
	{
		DimVector<3> tuv = _inverseMapping * Vect<3>(v1, domainPoint);
		RefPoint p(tuv(0), tuv(1), tuv(2));
		return p;
	}

	static GeometricMapping MappingInfo()
	{
		GeometricMapping mapping;
		mapping.NFunctions = 4; // 1, t, u, v

		mapping.Coeffs = vector<double>(mapping.NFunctions * mapping.NFunctions);
		for (int i = 0; i < mapping.NFunctions; i++)
			mapping.Coeffs[i * mapping.NFunctions + i] = 1; // Identity

		//                    t,u,v
		mapping.Exponents = { 0,0,0,   // 1   = t^0 * u^0 * v^0
							  1,0,0,   // t   = t^1 * u^0 * v^0
							  0,1,0,   // u   = t^0 * u^1 * v^0
							  0,0,1 }; // v   = t^0 * u^0 * v^1
		return mapping;
	}

	vector<double> MappingCoefficients() const override
	{
		// Refer to ConvertToDomain(). The basis functions are 1, t, u, v.
		//
		//          1,             t,            u
		return { v1.X, (v2.X - v1.X), (v3.X - v1.X), (v4.X - v1.X),   // X
				 v1.Y, (v2.Y - v1.Y), (v3.Y - v1.Y), (v4.Y - v1.Y),   // Y
				 v1.Z, (v2.Z - v1.Z), (v3.Z - v1.Z), (v4.Z - v1.Z) }; // Z
	}

	void Serialize(ostream& os) const override
	{
		os << "Tetrahedron";
		os << " ";
		v1.Serialize(os, 3);
		os << "--";
		v2.Serialize(os, 3);
		os << "--";
		v3.Serialize(os, 3);
		os << "--";
		v4.Serialize(os, 3);
	}

	//---------------------------------------------------------------------//
	// This is f***ing useless, it should be automatic due to inheritance! //
	// But without that it doesn't compile for some reason :-(             //
	//---------------------------------------------------------------------//

	double Integral(DomFunction globalFunction) const override
	{
		return PhysicalShapeWithConstantJacobian<3>::Integral(globalFunction);
	}
	double Integral(DomFunction globalFunction, int polynomialDegree) const override
	{
		return PhysicalShapeWithConstantJacobian<3>::Integral(globalFunction, polynomialDegree);
	}

	//-------------------------------------------------------------------//
	//                            Unit tests                             //
	//-------------------------------------------------------------------//

	static void Test()
	{
		DomPoint A(1, 1, 1);
		DomPoint B(4, 1, 1);
		DomPoint C(1, 4, 1);
		DomPoint D(1, 1, 4);

		Tetrahedron t(A, B, C, D);

		t.UnitTests();

		RefPoint ARef = t.ConvertToReference(A);
		assert(ARef == RefPoint(0, 0, 0));
		DomPoint ADom = t.ConvertToDomain(RefPoint(0, 0, 0));
		assert(A == ADom);

		RefPoint BRef = t.ConvertToReference(B);
		assert(BRef == RefPoint(1, 0, 0));
		DomPoint BDom = t.ConvertToDomain(RefPoint(1, 0, 0));
		assert(B == BDom);

		RefPoint CRef = t.ConvertToReference(C);
		assert(CRef == RefPoint(0, 1, 0));
		DomPoint CDom = t.ConvertToDomain(RefPoint(0, 1, 0));
		assert(C == CDom);

		RefPoint DRef = t.ConvertToReference(D);
		assert(DRef == RefPoint(0, 0, 1));
		DomPoint DDom = t.ConvertToDomain(RefPoint(0, 0, 1));
		assert(D == DDom);

		assert(abs(t.Measure() - pow(3, 3)/6.0) < 1e-14); // Tetra's volume is 1/6 of the cube's
	}
};

ReferenceTetrahedron Tetrahedron::RefTetra = ReferenceTetrahedron();