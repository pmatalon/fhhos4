#pragma once
#include <functional>
#include "../../Geometry/Point.h"
#include "Tensor.h"
#include "../../Mesh/PhysicalGroup.h"
using namespace std;

template <int Dim>
class DiffusionField
{
private:
	map<string, Tensor<Dim>> _tensorsByPhysicalPartName;
public:
	bool IsHomogeneous = true;
	bool IsIsotropic = true;

	DiffusionField()
	{}

	// This constructor is used when the domain is not partitioned into different physical parts.
	DiffusionField(const Tensor<Dim>& k)
	{
		this->_tensorsByPhysicalPartName.insert({ "domain", k });
		this->IsHomogeneous = true;
		this->IsIsotropic = k.IsIsotropic;
	}

	// This constructor is used when the domain is partitioned into 2 physical parts only.
	DiffusionField(string physicalPartName1, const Tensor<Dim>& k1, string physicalPartName2, const Tensor<Dim> k2)
	{
		this->_tensorsByPhysicalPartName.insert({ physicalPartName1, k1 });
		this->_tensorsByPhysicalPartName.insert({ physicalPartName2, k2 });
		this->IsHomogeneous = k1 == k2;
		this->IsIsotropic = k1.IsIsotropic && k2.IsIsotropic;
	}

	DiffusionField(const map<string, Tensor<Dim>>& tensors)
	{
		this->_tensorsByPhysicalPartName = tensors;
		this->IsHomogeneous = true;

		auto it = tensors.begin();
		const Tensor<Dim>& k1 = it->second;
		this->IsIsotropic = k1.IsIsotropic;
		it++;
		while (it != tensors.end())
		{
			const Tensor<Dim>& k = it->second;
			if (k != k1)
			{
				this->IsHomogeneous = false;
				if (!k.IsIsotropic)
					this->IsIsotropic = false;
			}
			it++;
		}
	}

	// Anisotropic field accross the whole domain.
	DiffusionField(double anisotropyRatio, double anisotropyAngle)
		: DiffusionField(Tensor<Dim>(1, anisotropyRatio, anisotropyAngle))
	{}


	Tensor<Dim>& ConstantDiffTensor(PhysicalGroup<Dim>* physicalPart)
	{
		if (physicalPart->Name.size() == 0)
			Utils::FatalError("A nameless physical part has been found. Please name all the physical parts of the geometry.");
		auto found = _tensorsByPhysicalPartName.find(physicalPart->Name);
		if (found == _tensorsByPhysicalPartName.end())
			Utils::FatalError("Unknown physical part '" + physicalPart->Name + "'");
		return found->second;
	}

	double LargestHeterogeneityRatio()
	{
		if (this->IsHomogeneous)
			return 1;

		double largestEigenValueMin = 1e30;
		double largestEigenValueMax = -1;
		for (auto it = _tensorsByPhysicalPartName.begin(); it != _tensorsByPhysicalPartName.end(); it++)
		{
			const Tensor<Dim>& k = it->second;
			largestEigenValueMin = min(largestEigenValueMin, k.LargestEigenValue);
			largestEigenValueMax = max(largestEigenValueMax, k.LargestEigenValue);
		}

		return largestEigenValueMax / largestEigenValueMin;
	}

	double LargestAnisotropyRatio()
	{
		double largestAnisotropyRatio = 1;
		for (auto it = _tensorsByPhysicalPartName.begin(); it != _tensorsByPhysicalPartName.end(); it++)
		{
			const Tensor<Dim>& k = it->second;
			largestAnisotropyRatio = max(largestAnisotropyRatio, k.AnisotropyRatio);
		}
		return largestAnisotropyRatio;
	}

	vector<const Tensor<Dim>*> Tensors()
	{
		vector<const Tensor<Dim>*> list;
		for (auto it = _tensorsByPhysicalPartName.begin(); it != _tensorsByPhysicalPartName.end(); it++)
			list.push_back(&it->second);
		return list;
	}
};