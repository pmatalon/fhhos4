#pragma once
#include <functional>
#include "../Geometry/Point.h"
#include "Tensor.h"
#include "../Mesh/PhysicalGroup.h"
using namespace std;

template <int Dim>
class DiffusionField
{
public:
	function<Tensor<Dim>*(PhysicalGroup<Dim>*)> ConstantDiffTensor = nullptr;

	bool IsHomogeneous = true;
	bool IsIsotropic = true;
	double HeterogeneityRatio = 1;

	Tensor<Dim>* K1 = nullptr;
	Tensor<Dim>* K2 = nullptr;

	DiffusionField()
	{}

	// This constructor is used when the domain is not partitioned into different physical parts.
	DiffusionField(Tensor<Dim>* k)
	{
		this->K1 = k;
		this->IsHomogeneous = true;
		this->IsIsotropic = k->IsIsotropic;
		this->HeterogeneityRatio = 1;
		this->ConstantDiffTensor = [k](PhysicalGroup<Dim>* physicalPart)
		{
			return k;
		};
	}

	// This constructor is used when the domain is partitioned into 2 physical parts only.
	DiffusionField(string physicalPartName1, Tensor<Dim>* k1, string physicalPartName2, Tensor<Dim>* k2)
	{
		this->K1 = k1;
		this->K2 = k2;
		this->IsHomogeneous = *K1 == *K2;
		this->IsIsotropic = K1->IsIsotropic && K2->IsIsotropic;
		this->HeterogeneityRatio = K1->LargestEigenValue / K2->LargestEigenValue;
		this->ConstantDiffTensor = [physicalPartName1, k1, physicalPartName2, k2](PhysicalGroup<Dim>* physicalPart)
		{
			if (physicalPart->Name.compare(physicalPartName1) == 0)
				return k1;
			else if (physicalPart->Name.compare(physicalPartName2) == 0)
				return k2;
			else if (physicalPart->Name.size() == 0)
				Utils::FatalError("A non-named physical part has been found. Please name all the physical parts of the geometry.");
			else
				Utils::FatalError("Unknown physical part '" + physicalPart->Name + "'");
		};
	}

	DiffusionField(const map<string, Tensor<Dim>*>& tensors)
	{
		this->IsHomogeneous = true;
		this->IsIsotropic = true;
		for (auto it = tensors.begin(); it != tensors.end(); it++)
		{
			Tensor<Dim>* k = it->second;
			if (!this->K1)
				this->K1 = k;
			else
			{
				bool sameTensors = *k == *this->K1;

				if (!this->K2 && !sameTensors)
					this->K2 = k;
				if (!sameTensors && k->LargestEigenValue != this->K1->LargestEigenValue)
					this->IsHomogeneous = false;
			}

			if (!k->IsIsotropic)
				this->IsIsotropic = false;
		}

		if (!this->IsHomogeneous)
		{
			if (K1->LargestEigenValue > K2->LargestEigenValue)
				this->HeterogeneityRatio = K1->LargestEigenValue / K2->LargestEigenValue;
			else
				this->HeterogeneityRatio = K2->LargestEigenValue / K1->LargestEigenValue;
		}

		this->ConstantDiffTensor = [tensors](PhysicalGroup<Dim>* physicalPart)
		{
			if (physicalPart->Name.size() == 0)
				Utils::FatalError("A non-named physical part has been found. Please name all the physical parts of the geometry.");
			auto found = tensors.find(physicalPart->Name);
			if (found == tensors.end())
				Utils::FatalError("Unknown physical part '" + physicalPart->Name + "'");
			return found->second;
		};
	}

	// Anisotropic field accross the whole domain.
	DiffusionField(double anisotropyRatio, double anisotropyAngle)
		: DiffusionField(new Tensor<Dim>(1, anisotropyRatio, anisotropyAngle))
	{}
};