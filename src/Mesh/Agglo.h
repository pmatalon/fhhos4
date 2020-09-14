#pragma once
#include "Element.h"
#include "Interface.h"
using namespace std;

template <int Dim>
class Agglo : public Element<Dim>
{
private:
	vector<Vertex*> _vertices;
	list<Element<Dim>*> _elemsToAgglomerate;
	vector<Face<Dim>*> _removedFaces;
public:
	Agglo(vector<Element<Dim>*> elements) :
		Element<Dim>(0),
		_vertices(elements[0]->Vertices()),
		_elemsToAgglomerate(elements.begin()+1, elements.end())
	{
		this->Faces = elements[0]->Faces;
		/*if (elements[0]->Number == 564)
		{
			cout << "%------------------------- Elements to agglomerate -------------------------------" << endl;
			for (Element<Dim>* e : elements)
			{
				cout << "% " << *e << endl;
				e->Shape()->ExportToMatlab();
			}
			cout << endl;
		}*/

		AddNeighboursRecursively(elements[0]);

		if (!_elemsToAgglomerate.empty())
		{
			cout << "%-------------------- Agglomerate -----------------" << endl;
			this->ExportToMatlab("b");
			for (Element<Dim>* e : _elemsToAgglomerate)
			{
				cout << "%-------------------- Element that hasn't been agglomerated -----------------" << endl;
				e->ExportToMatlab("r");
			}
			Utils::FatalError("Agglomeration failed: some elements have not been agglomerated for some reason...");
		}
	}

	inline vector<Vertex*> Vertices() const override
	{
		return _vertices;
	}

	inline vector<Face<Dim>*> RemovedFaces()
	{
		return _removedFaces;
	}

private:
	void AddNeighboursRecursively(Element<Dim>* e)
	{
		vector<Element<Dim>*> recursion;

		// Agglomerate the direct neighbours first
		for (Element<Dim>* neighbour : e->Neighbours())
		{
			typename list<Element<Dim>*>::iterator it;
			it = find(_elemsToAgglomerate.begin(), _elemsToAgglomerate.end(), neighbour);
			if (it != _elemsToAgglomerate.end()) // verify that the neighbour is in the list of elements to agglomerate
			{
				vector<Face<Dim>*> interfaceFaces = this->InterfaceWith(neighbour);
				Interface<Dim> interf(interfaceFaces);
				if (interf.HasHoles())
				{
					recursion.push_back(neighbour);
					continue; // we'll do this one later so we don't surround an element (one of its neighbours will process it)
				}

				_removedFaces = Utils::Join(_removedFaces, interfaceFaces);
				// Agglomeration
				/*if (e->Number == 534)
				{
					this->ExportToMatlab("b");
					cout << endl;
				}*/
				this->Agglomerate(neighbour, interfaceFaces);
				this->Faces = Utils::SymmetricDifference<Face<Dim>*>(this->Faces, neighbour->Faces);
				_elemsToAgglomerate.erase(it);
				recursion.push_back(neighbour);
			}
		}

		// And then agglomerate the neighbours of the neighbours
		for (Element<Dim>* neighbour : recursion)
			AddNeighboursRecursively(neighbour);
	}

	//-----------------------------------------------//
	// Dim-specific functions (implementation below) //
	//-----------------------------------------------//
private:
	void Agglomerate(Element<Dim>* e, const vector<Face<Dim>*>& interfaceFaces) { assert(false); }

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//
public:
	PhysicalShape<Dim>* Shape() const
	{
		assert(false);
	}

	virtual DimVector<Dim> OuterNormalVector(Face<Dim>* face) const override
	{
		assert(false);
	}

	void ExportToMatlab(string color = "r") const
	{
		MatlabScript script;
		script.PlotPolygonEdges(_vertices, color);
	}
};

template<>
void Agglo<2>::Agglomerate(Element<2>* e, const vector<Face<2>*>& interfaceFaces)
{
	/*if (e->Number == 534 && interfaceFaces.size() == 3)
	{
		MatlabScript script;
		cout << "%------------------------- Current agglomeration --------------------------------" << endl;
		script.PlotPolygonEdges(this->Vertices(), "m");
		cout << "%------------------------- Adding this element --------------------------------" << endl;
		script.PlotPolygonEdges(e->Vertices(), "c");
	}*/
	_vertices = PolygonalElement::MacroPolygonVertices(this, e, interfaceFaces);
	//cout << "%------------------------- Result --------------------------------" << endl;
	//script.PlotPolygonEdges(this->Vertices(), "r");
	//cout << endl;
}