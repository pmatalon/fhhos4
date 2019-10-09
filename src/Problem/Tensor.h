#pragma once
#include "../Utils/Types.h"
using namespace std;

template<int Dim>
class Tensor
{
public:
	Eigen::Matrix<double, Dim, Dim> TensorMatrix;
	double LargestEigenValue;
	double SmallestEigenValue;
	bool IsIsotropic;
	double AnisotropyRatio;
	double RotationAngle;

	Tensor(DimVector<Dim> anisotropyCoefficients)
		: Tensor(anisotropyCoefficients, 0)
	{ }

	Tensor(DimVector<Dim> anisotropyCoefficients, double rotationAngle)
	{
		//this->TensorMatrix = Utils::CreateTensor(anisotropyCoefficients, rotationAngle);
		Eigen::Matrix<double, Dim, Dim> Q;
		if (Dim == 1)
			Q << 1;
		else if (Dim == 2)
			Q << cos(rotationAngle), -sin(rotationAngle),
			sin(rotationAngle), cos(rotationAngle);
		else if (Dim == 3)
			Q << 1, 0, 0,
			0, cos(rotationAngle), -sin(rotationAngle),
			0, sin(rotationAngle), cos(rotationAngle);

		this->TensorMatrix = Q * anisotropyCoefficients.asDiagonal() * Q.transpose();

		this->LargestEigenValue = anisotropyCoefficients.maxCoeff();
		this->SmallestEigenValue = anisotropyCoefficients.minCoeff();
		this->IsIsotropic = LargestEigenValue == SmallestEigenValue;
		this->AnisotropyRatio = LargestEigenValue / SmallestEigenValue;
		this->RotationAngle = rotationAngle;
	}

	static Tensor<Dim>* Isotropic(double coeff)
	{
		DimVector<Dim> isotropyCoefficients = coeff * DimVector<Dim>::Ones(Dim);
		Tensor<Dim>* K = new Tensor<Dim>(isotropyCoefficients, 0);
		return K;
	}

	static Tensor<Dim>* Isotropic()
	{
		return Isotropic(1);
	}

	bool operator==(Tensor<Dim>& other)
	{
		return (this->TensorMatrix - other.TensorMatrix).isZero(0);
	}
};

template <int Dim>
DimVector<Dim> operator*(Tensor<Dim> const& K, DimVector<Dim> const& b)
{
	DimVector<Dim> result = K.TensorMatrix * b;
	return result;
}

template <int Dim>
DimVector<Dim> operator*(Tensor<Dim>* const K, DimVector<Dim> const& b)
{
	DimVector<Dim> result = K->TensorMatrix * b;
	return result;
}