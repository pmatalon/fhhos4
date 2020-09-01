#pragma once

enum class Action : unsigned
{
	None						= 0,
	ExtractSystem				= 1 << 1,
	ExtractComponentMatrices	= 1 << 2,
	ExportFaces					= 1 << 3,
	SolveSystem					= 1 << 4,
	ExportMultigridMatrices     = 1 << 5,
	ExtractSolution				= 1 << 6,
	LogAssembly					= 1 << 7,
	ComputeL2Error              = 1 << 8,
	UnitTests                   = 1 << 9,
	ExportSolutionToGMSH        = 1 << 10
};

Action operator &(Action lhs, Action rhs)
{
	return static_cast<Action> (
		static_cast<std::underlying_type<Action>::type>(lhs) &
		static_cast<std::underlying_type<Action>::type>(rhs)
		);
}

Action operator ^(Action lhs, Action rhs)
{
	return static_cast<Action> (
		static_cast<std::underlying_type<Action>::type>(lhs) ^
		static_cast<std::underlying_type<Action>::type>(rhs)
		);
}

Action operator ~(Action rhs)
{
	return static_cast<Action> (
		~static_cast<std::underlying_type<Action>::type>(rhs)
		);
}

Action& operator |=(Action &lhs, Action rhs)
{
	lhs = static_cast<Action> (
		static_cast<std::underlying_type<Action>::type>(lhs) |
		static_cast<std::underlying_type<Action>::type>(rhs)
		);

	return lhs;
}

Action& operator &=(Action &lhs, Action rhs)
{
	lhs = static_cast<Action> (
		static_cast<std::underlying_type<Action>::type>(lhs) &
		static_cast<std::underlying_type<Action>::type>(rhs)
		);

	return lhs;
}

Action& operator ^=(Action &lhs, Action rhs)
{
	lhs = static_cast<Action> (
		static_cast<std::underlying_type<Action>::type>(lhs)^
		static_cast<std::underlying_type<Action>::type>(rhs)
		);

	return lhs;
}