#pragma once
#include <fstream>
#include "Types.h"
#include "RotatingList.h"
#include "../Mesh/Vertex.h"
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

	void Add(string s)
	{
		Out() << s << endl;
	}

	void OpenFigure()
	{
		Out() << "figure, axes = gca; hold(axes, 'on');" << endl;
	}

	void PlotPoint(const DomPoint& p, string options = "k+")
	{
		Out() << "plot(axes, " << p.X << ", " << p.Y << ",'" << options << "');" << endl;
	}

	void Plot(const vector<Vertex*>& vertices, string options = "k+")
	{
		for (auto v : vertices)
			PlotPoint(*v, options);
	}

	void PlotSegment(Vertex* v1, Vertex* v2, string colorAndLineType, int lineWidth = 1)
	{
		PlotSegment(*v1, *v2, colorAndLineType, lineWidth);
	}

	void PlotSegment(const DomPoint& p1, const DomPoint& p2, string colorAndLineType, int lineWidth = 1)
	{
		Out() << "plot(axes, [" << p1.X << "; " << p2.X << "], [" << p1.Y << "; " << p2.Y << "], '" << colorAndLineType << "', 'LineWidth', " << lineWidth << ");" << endl;
		// use plot3 in 3D
	}

	void PlotPolygonEdges(const vector<Vertex*>& vertices, string colorAndLineType, int lineWidth = 1)
	{
		RotatingList<Vertex*> vert(vertices);
		for (int i = 0; i < vert.Size(); i++)
		{
			Vertex* v1 = vert.GetAndMoveNext();
			Vertex* v2 = vert.Get();
			PlotSegment(v1, v2, colorAndLineType, lineWidth);
		}
	}

	void PlotPolygonEdges(const vector<DomPoint>& vertices, string colorAndLineType, int lineWidth = 1)
	{
		RotatingList<DomPoint> vert(vertices);
		for (int i = 0; i < vert.Size(); i++)
		{
			DomPoint v1 = vert.GetAndMoveNext();
			DomPoint v2 = vert.Get();
			PlotSegment(v1, v2, colorAndLineType, lineWidth);
		}
	}

	void PlotTriangle(const DomPoint& A, const DomPoint& B, const DomPoint& C, string color)
	{
		Out() << "fill(axes, [" << A.X << ";" << B.X << ";" << C.X << "], [" << A.Y << "; " << B.Y << ";" << C.Y << "], '" << color << "', 'LineStyle','none', 'FaceAlpha', 0.5);" << endl;
		PlotSegment(A, B, color);
		PlotSegment(B, C, color);
		PlotSegment(C, A, color);
	}

	void PlotPolygon(const vector<DomPoint>& vertices, string color, string edgeLineType = "-")
	{
		stringstream X;
		stringstream Y;
		for (int i = 0; i < vertices.size(); i++)
		{
			X << vertices[i].X;
			Y << vertices[i].Y;
			if (i < vertices.size() - 1)
			{
				X << ";";
				Y << ";";
			}
		}
		Out() << "fill(axes, [" << X.str() << "], [" << Y.str() << "], '" << color << "', 'LineStyle','none', 'FaceAlpha', 0.5);" << endl;
		PlotPolygonEdges(vertices, color + edgeLineType);
	}

	void PlotText(const DomPoint& p, string text, string color = "k")
	{
		Out() << "text(axes, " << p.X << ", " << p.Y << ", '" << text << "', 'Color', '" << color << "');" << endl;
	}

	void Comment(string comment)
	{
		Out() << "% " << comment << endl;
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