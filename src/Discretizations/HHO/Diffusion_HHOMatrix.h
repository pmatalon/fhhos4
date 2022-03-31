#pragma once
#include "../../FunctionalBasis/FunctionalBasis.h"
#include "Diff_HHOElement.h"
#include "../../Utils/NonZeroCoefficients.h"

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

	inline BigNumber FirstRow(Diff_HHOElement<Dim>* e)
	{
		return e->Number() * _nDoFsPerElement;
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

	inline BigNumber FirstRow(Diff_HHOElement<Dim>* e)
	{
		return e->Number() * _nDoFsPerElement;
	}
	inline BigNumber FirstCol(Diff_HHOFace<Dim>* f)
	{
		return f->Number() * _nDoFsPerFace;
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

	inline BigNumber FirstRow(Diff_HHOFace<Dim>* f)
	{
		return f->Number() * _nDoFsPerFace;
	}

	inline BigNumber FirstCol(Diff_HHOFace<Dim>* f)
	{
		return FirstRow(f);
	}
};