#pragma once
#include <fstream>
#include "Types.h"
#include "RotatingList.h"
#include "../Geometry/Point.h"
using namespace std;

class MatlabScript
{
private:
	string _filePath;
	ofstream _file;

public:
	MatlabScript()
	{
		// sends to cout
	}
	MatlabScript(string filePath) : _filePath(filePath)
	{
		_file.open(filePath.c_str());
	}

	ostream& Out()
	{
		return _file.is_open() ? _file : std::cout;
	}

	void PlotPoint(DomPoint p, string options)
	{
		Out() << "plot(axes, " << p.X << ", " << p.Y << ",'" << options << "');" << endl;
	}

	void PlotSegment(Vertex* v1, Vertex* v2, string color)
	{
		PlotSegment(*v1, *v2, color);
	}

	void PlotSegment(DomPoint p1, DomPoint p2, string color)
	{
		Out() << "plot(axes, [" << p1.X << "; " << p2.X << "], [" << p1.Y << "; " << p2.Y << "], 'LineWidth', 1, 'Color', '" << color << "');" << endl;
		// use plot3 in 3D
	}

	void PlotPolygonEdges(const vector<Vertex*>& vertices, string color)
	{
		RotatingList<Vertex*> vert(vertices);
		for (int i = 0; i < vert.Size(); i++)
		{
			Vertex* v1 = vert.GetAndMoveNext();
			Vertex* v2 = vert.Get();
			PlotSegment(v1, v2, color);
		}
	}

	void PlotTriangle(DomPoint A, DomPoint B, DomPoint C, string color)
	{
		Out() << "fill(axes, [" << A.X << ";" << B.X << ";" << C.X << "], [" << A.Y << "; " << B.Y << ";" << C.Y << "], '" << color << "', 'LineStyle','none', 'FaceAlpha', 0.5);" << endl;
	}

	void PlotText(DomPoint p, string text, string color = "k")
	{
		Out() << "text(axes, " << p.X << ", " << p.Y << ", '" << text << "', 'Color', '" << color << "');" << endl;
	}

	void Close()
	{
		if (_file.is_open())
		{
			_file.close();
			cout << "Matlab script exported to \t" << _filePath << endl;
		}
	}

	~MatlabScript()
	{
		Close();
	}
};