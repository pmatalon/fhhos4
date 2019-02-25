#pragma once

enum class Action : unsigned
{
	None						= 0,
	ExtractSystem				= 1 << 1,
	ExtractComponentMatrices	= 1 << 2,
	ExtractMassMatrix			= 1 << 3,
	SolveSystem					= 1 << 4
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