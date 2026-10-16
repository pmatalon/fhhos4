#pragma once
#include "Smoother.h"
using namespace std;

class SmootherFactory
{
private:
	// The codes of the smoothers (-smoothers, see fhhos4 -h); the first one of a smoother is its name
	static const vector<pair<string, SmootherType>>& Codes()
	{
		static const vector<pair<string, SmootherType>> codes = {
			{ "j",      SmootherType::Jacobi },
			{ "j23",    SmootherType::DampedJacobi },
			{ "bj",     SmootherType::BlockJacobi },
			{ "bj23",   SmootherType::DampedBlockJacobi },
			{ "gs",     SmootherType::GaussSeidel },
			{ "rgs",    SmootherType::BackwardGaussSeidel },
			{ "sgs",    SmootherType::SymmetricGaussSeidel },
			{ "ags",    SmootherType::AlternatingGaussSeidel },
			{ "rags",   SmootherType::BackwardAlternatingGaussSeidel },
			{ "sor",    SmootherType::SOR },
			{ "rsor",   SmootherType::BackwardSOR },
			{ "bgs",    SmootherType::BlockGaussSeidel },
			{ "bsor",   SmootherType::BlockGaussSeidel },
			{ "rbgs",   SmootherType::BackwardBlockGaussSeidel },
			{ "rbsor",  SmootherType::BackwardBlockGaussSeidel },
			{ "sbgs",   SmootherType::SymmetricBlockGaussSeidel },
			{ "sbsor",  SmootherType::SymmetricBlockGaussSeidel },
			{ "abgs",   SmootherType::AlternatingBlockGaussSeidel },
			{ "absor",  SmootherType::AlternatingBlockGaussSeidel },
			{ "rabgs",  SmootherType::BackwardAlternatingBlockGaussSeidel },
			{ "rabsor", SmootherType::BackwardAlternatingBlockGaussSeidel },
			{ "hbgs",   SmootherType::HybridBlockGaussSeidel },
			{ "hrbgs",  SmootherType::HybridBackwardBlockGaussSeidel },
			{ "cheb",   SmootherType::ChebyshevFirstKind },
			{ "bcheb",  SmootherType::BlockChebyshevFirstKind },
			{ "cheb4",  SmootherType::ChebyshevFourthKind },
			{ "bcheb4", SmootherType::BlockChebyshevFourthKind }
		};
		return codes;
	}

	// Block SOR, or block Gauss-Seidel if omega = 1 (Gauss-Seidel if the blocks are scalar)
	static Smoother* BlockSOROrGaussSeidel(int blockSize, double omega, Direction direction, int nSmootherIterations)
	{
		if (blockSize == 1 && omega == 1)
			return new GaussSeidelSmoother(direction, nSmootherIterations);
		else if (omega == 1)
			return new BlockGaussSeidelSmoother(blockSize, direction, false, nSmootherIterations);
		else
			return new BlockSORSmoother(blockSize, omega, direction, nSmootherIterations);
	}

	static Smoother* HybridBlockGaussSeidel(int blockSize, Direction direction, int nSmootherIterations)
	{
		if (blockSize == 1) // the scalar kernels are faster
			return new GaussSeidelSmoother(direction, true, nSmootherIterations);
		return new BlockGaussSeidelSmoother(blockSize, direction, true, nSmootherIterations);
	}

public:
	// The smoother of the code; false if the code is unknown
	static bool FromCode(const string& code, SmootherType& smoother)
	{
		for (const auto& [c, s] : Codes())
		{
			if (c == code)
			{
				smoother = s;
				return true;
			}
		}
		return false;
	}

	// The code of the smoother (for the messages)
	static string Code(SmootherType smoother)
	{
		if (smoother == SmootherType::Default)
			return "<symmetric smoother>";
		for (const auto& [c, s] : Codes())
		{
			if (s == smoother)
				return c;
		}
		return "<unknown smoother " + to_string((unsigned)smoother) + ">";
	}

	// The Chebyshev smoothers: symmetric, for a symmetric positive definite matrix
	static bool IsChebyshev(SmootherType smoother)
	{
		return smoother == SmootherType::ChebyshevFirstKind || smoother == SmootherType::BlockChebyshevFirstKind ||
		       smoother == SmootherType::ChebyshevFourthKind || smoother == SmootherType::BlockChebyshevFourthKind;
	}

	static Smoother* Create(SmootherType smoother, int nSmootherIterations, int blockSize, double omega)
	{
		switch (smoother)
		{
			case SmootherType::GaussSeidel:
				return new GaussSeidelSmoother(Direction::Forward, nSmootherIterations);
			case SmootherType::BackwardGaussSeidel:
				return new GaussSeidelSmoother(Direction::Backward, nSmootherIterations);
			case SmootherType::SymmetricGaussSeidel:
				return new GaussSeidelSmoother(Direction::Symmetric, nSmootherIterations);
			case SmootherType::AlternatingGaussSeidel:
				return new GaussSeidelSmoother(Direction::AlternatingForwardFirst, nSmootherIterations);
			case SmootherType::BackwardAlternatingGaussSeidel:
				return new GaussSeidelSmoother(Direction::AlternatingBackwardFirst, nSmootherIterations);

			case SmootherType::SOR:
				return BlockSOROrGaussSeidel(1, omega, Direction::Forward, nSmootherIterations);
			case SmootherType::BackwardSOR:
				return BlockSOROrGaussSeidel(1, omega, Direction::Backward, nSmootherIterations);

			case SmootherType::BlockGaussSeidel:
				return BlockSOROrGaussSeidel(blockSize, omega, Direction::Forward, nSmootherIterations);
			case SmootherType::BackwardBlockGaussSeidel:
				return BlockSOROrGaussSeidel(blockSize, omega, Direction::Backward, nSmootherIterations);
			case SmootherType::SymmetricBlockGaussSeidel:
				return BlockSOROrGaussSeidel(blockSize, omega, Direction::Symmetric, nSmootherIterations);
			case SmootherType::AlternatingBlockGaussSeidel:
				return BlockSOROrGaussSeidel(blockSize, omega, Direction::AlternatingForwardFirst, nSmootherIterations);
			case SmootherType::BackwardAlternatingBlockGaussSeidel:
				return BlockSOROrGaussSeidel(blockSize, omega, Direction::AlternatingBackwardFirst, nSmootherIterations);

			case SmootherType::HybridBlockGaussSeidel:
				return HybridBlockGaussSeidel(blockSize, Direction::Forward, nSmootherIterations);
			case SmootherType::HybridBackwardBlockGaussSeidel:
				return HybridBlockGaussSeidel(blockSize, Direction::Backward, nSmootherIterations);

			case SmootherType::Jacobi:
				return new BlockJacobiSmoother(1, omega, nSmootherIterations);
			case SmootherType::DampedJacobi:
				return new BlockJacobiSmoother(1, 2.0/3.0, nSmootherIterations);
			case SmootherType::BlockJacobi:
				return new BlockJacobiSmoother(blockSize, omega, nSmootherIterations);
			case SmootherType::DampedBlockJacobi:
				return new BlockJacobiSmoother(blockSize, 2.0/3.0, nSmootherIterations);

			case SmootherType::ChebyshevFirstKind:
				return new ChebyshevSmoother(1, ChebyshevKind::First, nSmootherIterations);
			case SmootherType::BlockChebyshevFirstKind:
				return new ChebyshevSmoother(blockSize, ChebyshevKind::First, nSmootherIterations);
			case SmootherType::ChebyshevFourthKind:
				return new ChebyshevSmoother(1, ChebyshevKind::Fourth, nSmootherIterations);
			case SmootherType::BlockChebyshevFourthKind:
				return new ChebyshevSmoother(blockSize, ChebyshevKind::Fourth, nSmootherIterations);

			case SmootherType::Default: // resolved before (AlgebraicSolverFactory::SetMultigridParameters())
				break;
		}
		Utils::FatalError("Unmanaged smoother: '" + Code(smoother) + "'");
		return nullptr;
	}
};
