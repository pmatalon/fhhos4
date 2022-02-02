#pragma once
#include "../../Utils/Types.h"
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

	// Isotropic tensor
	Tensor()
		: Tensor(1)
	{ }

	// Isotropic tensor with kappa as diffusion coefficient
	Tensor(double kappa)
		: Tensor(kappa * DimVector<Dim>::Ones(), 0)
	{ }

	// Anisotropic tensor with no rotation angle
	Tensor(DimVector<Dim> anisotropyCoefficients)
		: Tensor(anisotropyCoefficients, 0)
	{ }

	Tensor(double kappa, double anisotropyRatio, double rotationAngle)
	{
		DimVector<Dim> anisotropyCoefficients = kappa * DimVector<Dim>::Ones();
		anisotropyCoefficients[0] = anisotropyRatio * anisotropyCoefficients[0];
		Init(anisotropyCoefficients, rotationAngle);
	}

	// Anisotropic tensor.
	// 'rotationAngle' must be expressed in radians.
	Tensor(DimVector<Dim> anisotropyCoefficients, double rotationAngle)
	{
		Init(anisotropyCoefficients, rotationAngle);
	}

private:
	void Init(DimVector<Dim> anisotropyCoefficients, double rotationAngle)
	{
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

public:
	friend bool operator==(const Tensor<Dim>& k1, const Tensor<Dim>& k2)
	{
		return (k1.TensorMatrix - k2.TensorMatrix).isZero(0);
	}

	friend bool operator!=(const Tensor<Dim>& k1, const Tensor<Dim>& k2)
	{
		return !(k1 == k2);
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