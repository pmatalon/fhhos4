#pragma once
#include <functional>
#include "../Geometry/Point.h"
#include "Tensor.h"
using namespace std;

template <int Dim>
class DiffusionField
{
public:
	string Partition;
	function<bool(DomPoint)> IsInPart1 = nullptr;
	Tensor<Dim>* K1;
	Tensor<Dim>* K2;
	bool IsHomogeneous = true;
	bool IsIsotropic = true;
	double HeterogeneityRatio = 1;

	DiffusionField(string partition, Tensor<Dim>* k1, Tensor<Dim>* k2)
	{
		this->Partition = partition;
		if (partition.compare("halves") == 0 || Dim == 1)
			this->IsInPart1 = [](DomPoint p) { return p.X < 0.5; };
		else if (partition.compare("chiasmus") == 0)
			this->IsInPart1 = [](DomPoint p) { return (p.X < 0.5 && p.Y >= 0.5) || (p.X >= 0.5 && p.Y < 0.5); };
		else if (partition.compare("circle") == 0)
			this->IsInPart1 = [](DomPoint p) { return pow(p.X - 0.5, 2) + pow(p.Y - 0.5, 2) <= pow(0.25, 2); };
		else
			assert(false);

		this->K1 = k1;
		this->K2 = k2;
		this->IsHomogeneous = *K1 == *K2;
		this->IsIsotropic = K1->IsIsotropic && K2->IsIsotropic;
		this->HeterogeneityRatio = K1->LargestEigenValue / K2->LargestEigenValue;
	}

	Tensor<Dim>* DiffTensor(DomPoint p)
	{
		return this->IsInPart1(p) ? this->K1 : this->K2;
	}
};