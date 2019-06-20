#pragma once
#include <algorithm>
#include "../CartesianElement.h"
#include "Edge.h"
using namespace std;

class RectangularPolygon : public CartesianElement<2>
{
public:
	vector<Face<2>*> NorthFaces;
	vector<Face<2>*> SouthFaces;
	vector<Face<2>*> EastFaces;
	vector<Face<2>*> WestFaces;

	DomPoint BottomLeftCorner;
	DomPoint TopLeftCorner;
	DomPoint TopRightCorner;
	DomPoint BottomRightCorner;

	RectangularPolygon(int number, DomPoint bottomLeftCorner, double widthX, double widthY) :
		Element(number), 
		CartesianElement(number, bottomLeftCorner, widthX, widthY),
		BottomLeftCorner(bottomLeftCorner),
		TopLeftCorner(bottomLeftCorner.X, bottomLeftCorner.Y + widthY),
		TopRightCorner(bottomLeftCorner.X + widthX, bottomLeftCorner.Y + widthY),
		BottomRightCorner(bottomLeftCorner.X + widthX, bottomLeftCorner.Y)
	{}

	void AddNorthFace(Face<2>* face)
	{
		this->AddFace(face);
		this->NorthFaces.push_back(face);
	}

	void AddSouthFace(Face<2>* face)
	{
		this->AddFace(face);
		this->SouthFaces.push_back(face);
	}

	void AddEastFace(Face<2>* face)
	{
		this->AddFace(face);
		this->EastFaces.push_back(face);
	}

	void AddWestFace(Face<2>* face)
	{
		this->AddFace(face);
		this->WestFaces.push_back(face);
	}

	void SetWestFacesFromNeighbour(RectangularPolygon* leftNeighbour)
	{
		for (Face<2>* f : leftNeighbour->EastFaces)
		{
			AddWestFace(f);
			f->Element2 = this;
		}
	}

	void SetSouthFacesFromNeighbour(RectangularPolygon* bottomNeighbour)
	{
		for (Face<2>* f : bottomNeighbour->NorthFaces)
		{
			AddSouthFace(f);
			f->Element2 = this;
		}
	}

	//-------------------------------------------------------//
	//                 Element implementation                //
	//-------------------------------------------------------//

	vector<double> OuterNormalVector(Face<2>* face)
	{
		if (isIn(this->NorthFaces, face))
			return vector<double>{ 0, 1 };
		if (isIn(this->SouthFaces, face))
			return vector<double>{ 0, -1 };
		if (isIn(this->EastFaces, face))
			return vector<double>{ 1, 0 };
		if (isIn(this->WestFaces, face))
			return vector<double>{ -1, 0 };
		assert(false);
	}

private:
	bool isIn(vector<Face<2>*> faces, Face<2>* f)
	{
		auto it = find(faces.begin(), faces.end(), f);
		if (it != faces.end())
			return true;
		return false;
	}

};