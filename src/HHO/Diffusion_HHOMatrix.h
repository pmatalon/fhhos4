#pragma once
#include "../FunctionalBasis/FunctionalBasis.h"
#include "Diff_HHOElement.h"
#include "../Utils/NonZeroCoefficients.h"

template <int Dim>
class A_T_T_Block : public NonZeroCoefficients
{
private:
	int _nDoFsPerElement;
public:
	A_T_T_Block() {}
	A_T_T_Block(int nDoFsPerElement)
	{
		_nDoFsPerElement = nDoFsPerElement;
	}

	inline BigNumber FirstRow(Element<Dim>* e)
	{
		return e->Number * _nDoFsPerElement;
	}

	inline BigNumber Row(Element<Dim>* e, BasisFunction<Dim>* cellPhi)
	{
		return FirstRow(e) + cellPhi->LocalNumber;
	}
	inline BigNumber Col(Element<Dim>* e, BasisFunction<Dim>* cellPhi)
	{
		return Row(e, cellPhi);
	}
};




template <int Dim>
class A_T_F_Block : public NonZeroCoefficients
{
private:
	int _nDoFsPerElement;
	int _nDoFsPerFace;
public:
	A_T_F_Block() {}
	A_T_F_Block(int nDoFsPerElement, int nDoFsPerFace)
	{
		_nDoFsPerElement = nDoFsPerElement;
		_nDoFsPerFace = nDoFsPerFace;
	}

	inline BigNumber FirstRow(Element<Dim>* e)
	{
		return e->Number * _nDoFsPerElement;
	}
	inline BigNumber FirstCol(Face<Dim>* f)
	{
		return f->Number * _nDoFsPerFace;
	}

	inline BigNumber Row(Element<Dim>* e, BasisFunction<Dim>* cellPhi)
	{
		return FirstRow(e) + cellPhi->LocalNumber;
	}
	inline BigNumber Col(Face<Dim>* f, BasisFunction<Dim-1>* facePhi)
	{
		return FirstCol(f) + facePhi->LocalNumber;
	}
};



template <int Dim>
class A_F_F_Block : public NonZeroCoefficients
{
private:
	int _nDoFsPerFace;
public:
	A_F_F_Block() {}
	A_F_F_Block(int nDoFsPerFace)
	{
		_nDoFsPerFace = nDoFsPerFace;
	}

	inline BigNumber FirstRow(Face<Dim>* f)
	{
		return f->Number * _nDoFsPerFace;
	}

	inline BigNumber Row(Face<Dim>* f, BasisFunction<Dim-1>* facePhi)
	{
		return FirstRow(f) + facePhi->LocalNumber;
	}
	inline BigNumber Col(Face<Dim>* f, BasisFunction<Dim-1>* facePhi)
	{
		return Row(f, facePhi);
	}
};