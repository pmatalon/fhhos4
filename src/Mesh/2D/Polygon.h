#pragma once
#include "../../DG/Diff_DGElement.h"
#include "../../HHO/Diff_HHOElement.h"
#include "PolygonalShape.h"
#include "../../Utils/RotatingList.h"
using namespace std;

class Polygon : public Diff_DGElement<2>, public Diff_HHOElement<2>
{
private:
	PolygonalShape* _shape = nullptr;

public:
	// Constructor creating the polygon from the adjonction of two elements
	Polygon(int number, Element<2>* e1, Element<2>* e2, vector<Face<2>*> facesToRemove) :
		Element(number),
		Diff_DGElement<2>(number),
		Diff_HHOElement<2>(number)
	{
		_shape = new PolygonalShape(MacroPolygonVertices(e1, e2, facesToRemove));
	}

	Polygon(int number, vector<Vertex*> vertices) :
		Element(number),
		Diff_DGElement<2>(number),
		Diff_HHOElement<2>(number)
	{
		_shape = new PolygonalShape(vertices);
	}

	inline vector<Vertex*> Vertices()
	{
		return _shape->Vertices();
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	GeometricShapeWithReferenceShape<2>* Shape() const
	{
		return _shape;
	}

	DimVector<2> OuterNormalVector(Face<2>* face) const
	{
		DimVector<2> n;
		Vertex* A = face->Vertices()[0];
		Vertex* B = face->Vertices()[1];

		// Condition 1: n.AB = 0
		// =>  n = (-AB.Y, AB.X)
		n << A->Y - B->Y, B->X - A->X;

		// Condition 2: n.AC < 0
		DimVector<2> AC = this->Center() - *A;
		if (n.dot(AC) > 0)
			n = -1 * n;
		
		n = n.normalized();
		return n;
	}

	static vector<Vertex*> MacroPolygonVertices(Element<2>* e1, Element<2>* e2, vector<Face<2>*> facesToRemove)
	{
		// We need to keep them in the direct order!
		vector<Vertex*> macroElementVertices;

		if (!e1)
		{
			assert(e2);
			macroElementVertices = e2->Vertices();
		}
		else if (!e2)
		{
			assert(e1);
			macroElementVertices = e1->Vertices();
		}
		else
		{
			assert(!facesToRemove.empty());

			// Find the first and last vertices (in the direct order) of the interface between the two elements
			RotatingList<Vertex*> e1Vertices(e1->Vertices());
			Vertex* firstInterfaceVertex = nullptr;
			Vertex* lastInterfaceVertex = nullptr;
			int iFirstInterfaceVertex = -1;
			int iLastInterfaceVertex = -1;
			int n = Face<2>::NumberOfFacesContainingVertex(facesToRemove, e1Vertices.Get());
			while (n == 1)
			{
				e1Vertices.MoveNext();
				n = Face<2>::NumberOfFacesContainingVertex(facesToRemove, e1Vertices.Get());
			}
			if (n == 0)
			{
				e1Vertices.MoveNext();
				while (!Face<2>::IsInFaces(facesToRemove, e1Vertices.Get()))
					e1Vertices.MoveNext();
				iFirstInterfaceVertex = e1Vertices.Index();
				firstInterfaceVertex = e1Vertices.Get();

				e1Vertices.MoveNext();
				while (Face<2>::IsInTwoFaces(facesToRemove, e1Vertices.Get()))
					e1Vertices.MoveNext();
				iLastInterfaceVertex = e1Vertices.Index();
				lastInterfaceVertex = e1Vertices.Get();
			}
			else // n >= 2
			{
				e1Vertices.MoveNext();
				while (Face<2>::IsInTwoFaces(facesToRemove, e1Vertices.Get()))
					e1Vertices.MoveNext();
				iLastInterfaceVertex = e1Vertices.Index();
				lastInterfaceVertex = e1Vertices.Get();

				e1Vertices.MoveNext();
				while (!Face<2>::IsInFaces(facesToRemove, e1Vertices.Get()))
					e1Vertices.MoveNext();
				iFirstInterfaceVertex = e1Vertices.Index();
				firstInterfaceVertex = e1Vertices.Get();
			}


			// Add all the vertices of e1 from lastInterfaceVertex to firstInterfaceVertex
			e1Vertices.GoTo(iLastInterfaceVertex);
			while (e1Vertices.Index() != iFirstInterfaceVertex)
			{
				macroElementVertices.push_back(e1Vertices.Get());
				e1Vertices.MoveNext();
			}
			macroElementVertices.push_back(firstInterfaceVertex);


			// Skip all vertices of e2 until you get to firstInterfaceVertex
			RotatingList<Vertex*> e2Vertices(e2->Vertices());
			e2Vertices.GoTo(firstInterfaceVertex);

			// Add e2's vertices
			e2Vertices.MoveNext();
			while (e2Vertices.Get() != lastInterfaceVertex)
			{
				macroElementVertices.push_back(e2Vertices.Get());
				e2Vertices.MoveNext();
			}
		}

		return macroElementVertices;
	}

	virtual ~Polygon()
	{
		if (_shape)
			delete _shape;
	}

};