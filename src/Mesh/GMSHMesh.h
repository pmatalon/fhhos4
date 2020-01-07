#pragma once
#include <gmsh.h>
#include "PolyhedralMesh.h"
#include "2D/Quadrilateral.h"
using namespace std;

enum GMSHElementTypes
{
	GMSH_Segment = 1,
	GMSH_Triangle = 2,
	GMSH_Quadrilateral = 3,
	GMSH_Tetrahedron = 4
};

template <int Dim>
class MeshVertex : public Vertex
{
public:
	vector<Element<Dim>*> Elements;
	vector<Face<Dim>*> Faces;

	MeshVertex(BigNumber number, double x) : Vertex(number, x) {}
	MeshVertex(BigNumber number, double x, double y) : Vertex(number, x, y) {}
	MeshVertex(BigNumber number, double x, double y, double z) : Vertex(number, x, y, z) {}
};

template <int Dim>
class GMSHMesh : public PolyhedralMesh<Dim>
{
private:
	map<size_t, BigNumber> _elementExternalNumbers;
	map<size_t, BigNumber> _vertexExternalNumbers;
public:
	GMSHMesh(string mshFile) : PolyhedralMesh<Dim>()
	{
		gmsh::initialize();
		gmsh::open(mshFile);

		for (int i=0; i<5; i++)
			gmsh::model::mesh::refine();

		//----------//
		// Vertices //
		//----------//

		vector<size_t> nodeTags;
		vector<double> coord;
		vector<double> parametricCoord;
		bool returnParametricCoord = false;
		bool includeBoundary = true;
		gmsh::model::mesh::getNodes(nodeTags, coord, parametricCoord, Dim, -1, includeBoundary, returnParametricCoord);

		this->Vertices.reserve(nodeTags.size());
		BigNumber vertexNumber = 0;
		for (size_t i = 0; i < nodeTags.size(); i++)
		{
			if (_vertexExternalNumbers.find(nodeTags[i]) == _vertexExternalNumbers.end())
			{
				MeshVertex<Dim>* v = new MeshVertex<Dim>(vertexNumber, coord[3 * i], coord[3 * i + 1], coord[3 * i + 2]);
				_vertexExternalNumbers.insert({ nodeTags[i], vertexNumber });
				this->Vertices.push_back(v);
				vertexNumber++;
			}
		}

		//----------//
		// Elements //
		//----------//

		vector<int> elementTypes;
		vector<vector<size_t>> elementTags;
		vector<vector<size_t>> elemNodeTags;
		gmsh::model::mesh::getElements(elementTypes, elementTags, elemNodeTags, Dim, -1);

		BigNumber elemNumber = 0;
		for (int i = 0; i < elementTypes.size(); i++)
		{
			int elemType = elementTypes[i];
			vector<size_t> elements = elementTags[i];
			vector<size_t> elementNodes = elemNodeTags[i];

			size_t elemNodeIndex = 0;
			for (size_t j = 0; j < elements.size(); j++)
			{
				Element<Dim>* e = nullptr;
				if (elemType == GMSH_Quadrilateral)
				{
					size_t nodeTag1 = elementNodes[elemNodeIndex];
					size_t nodeTag2 = elementNodes[elemNodeIndex + 1];
					size_t nodeTag3 = elementNodes[elemNodeIndex + 2];
					size_t nodeTag4 = elementNodes[elemNodeIndex + 3];

					e = new Quadrilateral(elemNumber, GetVertex(nodeTag1), GetVertex(nodeTag2), GetVertex(nodeTag3), GetVertex(nodeTag4));
					_elementExternalNumbers.insert({ elements[j], elemNumber });
					elemNumber++;
					elemNodeIndex += 4;
				}
				else if (elemType == GMSH_Triangle)
				{
					size_t nodeTag1 = elementNodes[elemNodeIndex];
					size_t nodeTag2 = elementNodes[elemNodeIndex + 1];
					size_t nodeTag3 = elementNodes[elemNodeIndex + 2];

					e = new Triangle(elemNumber, GetVertex(nodeTag1), GetVertex(nodeTag2), GetVertex(nodeTag3));
					_elementExternalNumbers.insert({ elements[j], elemNumber });
					elemNumber++;
					elemNodeIndex += 3;
				}
				else
					assert(false && "GMSH element type not managed.");

				for (Vertex* v : e->Shape()->Vertices())
				{
					MeshVertex<Dim>* mv = (MeshVertex<Dim>*)v;
					mv->Elements.push_back(e);
				}
				this->Elements.push_back(e);
			}
		}

		//-----------//
		//   Faces   //
		//-----------//

		BigNumber faceNumber = 0;
		for (int i = 0; i < elementTypes.size(); i++)
		{
			int elemType = elementTypes[i];

			if (Dim == 2)
			{
				vector<size_t> edgeNodes;
				bool onlyPrimaryNodes = true;
				gmsh::model::mesh::getElementEdgeNodes(elemType, edgeNodes, -1, onlyPrimaryNodes);

				for (size_t j = 0; j < edgeNodes.size(); j+=2)
				{
					MeshVertex<Dim>* v1 = (MeshVertex<Dim>*)GetVertex(edgeNodes[j]);
					MeshVertex<Dim>* v2 = (MeshVertex<Dim>*)GetVertex(edgeNodes[j + 1]);

					bool edgeAlreadyExists = false;
					for (Face<Dim>* f : v1->Faces)
					{
						Edge* edge = dynamic_cast<Edge*>(f);
						if (edge->Vertex1() == v2 || edge->Vertex2() == v2)
						{
							edgeAlreadyExists = true;
							break;
						}
					}
					if (edgeAlreadyExists)
						continue;

					vector<Element<Dim>*> neighbours;
					for (Element<Dim>* e1 : v1->Elements)
					{
						for (Element<Dim>* e2 : v2->Elements)
						{
							if (e1 == e2)
							{
								neighbours.push_back(e1);
								break;
							}
						}
					}

					Edge* edge;
					if (neighbours.size() == 1)
					{
						edge = new Edge(faceNumber++, v1, v2, neighbours[0]);
						neighbours[0]->AddFace(edge);
						this->BoundaryFaces.push_back(edge);
					}
					else if (neighbours.size() == 2)
					{
						edge = new Edge(faceNumber++, v1, v2, neighbours[0], neighbours[1]);
						neighbours[0]->AddFace(edge);
						neighbours[1]->AddFace(edge);
						this->InteriorFaces.push_back(edge);
					}
					else
						assert(false);
					this->Faces.push_back(edge);
					v1->Faces.push_back(edge);
					v2->Faces.push_back(edge);
				}
			}
		}
	}

private:
	inline Vertex* GetVertex(BigNumber nodeTag)
	{
		BigNumber internalNumber = _vertexExternalNumbers.at(nodeTag);
		return this->Vertices[internalNumber];
	}

public:
	~GMSHMesh()
	{
		gmsh::finalize();
	}

	string Description() override
	{
		return "GMSH";
	}

	string FileNamePart() override
	{
		return "GMSH";
	}

	double H() override
	{
		
	}

	void CoarsenMesh(CoarseningStrategy strategy) override
	{
		assert(false);
	}

};