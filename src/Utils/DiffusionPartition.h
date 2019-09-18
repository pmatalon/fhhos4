#pragma once
#include <functional>
#include "../Mesh/Point.h"
#include "Tensor.h"
using namespace std;

template <int Dim>
class DiffusionPartition
{
public:
	string Partition;
	function<bool(DomPoint)> IsInPart1 = nullptr;
	double Kappa1 = 1; // still used in DG
	double Kappa2 = 1; // still used in DG
	Tensor<Dim>* K1;
	Tensor<Dim>* K2;
	bool IsHomogeneous = true;
	bool IsIsotropic = true;
	double HeterogeneityRatio = 1;

	DiffusionPartition(string partition, Tensor<Dim>* k1, Tensor<Dim>* k2)
	{
		this->Partition = partition;
		if (partition.compare("halves") == 0 || Dim == 1)
			this->IsInPart1 = [](DomPoint p) { return p.X < 0.5; };
		else if (partition.compare("chiasmus") == 0)
			this->IsInPart1 = [](DomPoint p) { return (p.X < 0.5 && p.Y >= 0.5) || (p.X >= 0.5 && p.Y < 0.5); };
		else
			assert(false);

		this->K1 = k1;
		this->K2 = k2;
		this->Kappa1 = k1->LargestEigenValue; // still used in DG
		this->Kappa2 = k2->LargestEigenValue; // still used in DG
		this->IsHomogeneous = *K1 == *K2;
		this->IsIsotropic = K1->IsIsotropic && K2->IsIsotropic;
		this->HeterogeneityRatio = K1->LargestEigenValue / K2->LargestEigenValue;
	}

	// Still used in DG
	double Coefficient(DomPoint p)
	{
		return this->IsInPart1(p) ? this->Kappa1 : this->Kappa2;
	}

	Tensor<Dim>* DiffTensor(DomPoint p)
	{
		return this->IsInPart1(p) ? this->K1 : this->K2;
	}
};