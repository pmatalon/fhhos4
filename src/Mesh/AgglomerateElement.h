#pragma once
#include "../DG/Diff_DGElement.h"
#include "../HHO/Diff_HHOElement.h"
#include "AgglomerateShape.h"
using namespace std;

template <int Dim>
class AgglomerateElement : public Diff_DGElement<Dim>, public Diff_HHOElement<Dim>
{
private:
	AgglomerateShape<Dim>* _shape;
public:
	AgglomerateElement(int number, Element<Dim>* e1) :
		Element<Dim>(number),
		Diff_DGElement<Dim>(number),
		Diff_HHOElement<Dim>(number)
	{
		_shape = new AgglomerateShape<Dim>(e1->Shape());
	}

	void Add(Element<Dim>* e)
	{
		_shape->Add(e->Shape());
	}

	void Init()
	{
		assert(this->Faces.size() > 2);

		set<Vertex*> vertices;
		for (Face<Dim>* f : this->Faces)
		{
			for (Vertex* v : f->Vertices())
				vertices.insert(v);
		}
		
		_shape->Init(vector<Vertex*>(vertices.begin(), vertices.end()));
	}

	inline vector<Vertex*> Vertices()
	{
		return _shape->Vertices();
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	GeometricShapeWithReferenceShape<Dim>* Shape() const
	{
		return _shape;
	}

	virtual DimVector<Dim> OuterNormalVector(Face<Dim>* face) const override
	{
		assert(false);
	}

private:
	virtual void RemoveIntersections(vector<Face<Dim>*> faces, Face<Dim>* mergedFace) override
	{
		assert(false);
	}

public:
	virtual ~AgglomerateElement()
	{
		if (_shape)
			delete _shape;
	}

};

template<>
DimVector<2> AgglomerateElement<2>::OuterNormalVector(Face<2>* face) const
{
	DimVector<2> n;
	Vertex* A = face->Vertices()[0];
	Vertex* B = face->Vertices()[1];

	// Condition 1: n.AB = 0
	// =>  n = (-AB.Y, AB.X)
	n << A->Y - B->Y, B->X - A->X;

	n = n.normalized();

	// Condition 2: n.AC < 0
	// Find a 3rd point C inside the element to implement n.AC < 0
	// (the element's center doesn't work because the polygon might not be convex)
	GeometricShapeWithReferenceShape<2>* ss = _shape->ClosestSubShape(face->Center());
	DomPoint C = ss->Center();
	DimVector<2> AC = Vect<2>(A, C);
	if (n.dot(AC) > 0)
		n = -1 * n;

	return n;
}

template<>
DimVector<3> AgglomerateElement<3>::OuterNormalVector(Face<3>* face) const
{
	assert(false && "To be implemented");
}

template<>
void AgglomerateElement<2>::RemoveIntersections(vector<Face<2>*> oldFaces, Face<2>* newFace)
{
	// Identification of the vertices to remove (those which belong to 2 faces)
	vector<Vertex*> verticesToRemove;
	for (int i = 0; i < oldFaces.size(); i++)
	{
		Face<2>* fi = oldFaces[i];
		for (Vertex* v : fi->Vertices())
		{
			for (int j = i + 1; j < oldFaces.size(); j++)
			{
				Face<2>* fj = oldFaces[j];
				if (fj->HasVertex(v))
				{
					verticesToRemove.push_back(v);
					break;
				}
			}
		}
	}

	// Adds all the vertices in the same order, except the intersection vertices
	vector<Vertex*> newVertices;
	for (Vertex* v : _shape->Vertices())
	{
		if (!v->IsIn(verticesToRemove))
			newVertices.push_back(v);
	}

	assert(newVertices.size() > 2);

	// Reshape the shape:
	// Create a new virtual vertex at the center of the face, and move all the vertices located in the interior of the old faces to that point
	Vertex* C = new Vertex(0, newFace->Center());

	set<Vertex*> verticesToMove;
	for (auto ss : _shape->SubShapes())
	{
		for (auto v : ss->Vertices())
		{
			if (newFace->HasVertex(v))
				continue;
			for (auto f : oldFaces)
			{
				if (f->HasVertex(v) || f->Contains(*v))
				{
					verticesToMove.insert(v);
					break;
				}
			}
		}
	}

	for (Vertex* oldV : verticesToMove)
		_shape->ReshapeByMovingIntersection(oldV, C);

	_shape->CheckIfTwoSubShapesOverlap();

	// Adds all the vertices in the same order, except the intersection vertices, and add the center of the new face
	/*vector<Vertex*> newVertices;
	bool newFaceCenterAdded = false;
	for (Vertex* v : _shape->Vertices())
	{
		if (!v->IsIn(verticesToMove))
			newVertices.push_back(v);
		else if (!newFaceCenterAdded)
		{
			newVertices.push_back(C);
			newFaceCenterAdded = true;
		}
	}*/

	//assert(newVertices.size() > 2);

	// Recompute all the information about the shape (diameter, etc.)
	_shape->Init(newVertices);
}

template<>
void AgglomerateElement<3>::RemoveIntersections(vector<Face<3>*> faces, Face<3>* newFace)
{
	assert(false && "To be implemented");
}