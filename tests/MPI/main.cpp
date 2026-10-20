// The tests on several processes: run by ctest with mpiexec on 1 to 4 processes (MPI_np1...4). The results printed by
// the first process; a test fails if it fails on one process.
#include <gtest/gtest.h>
#include "Utils/Communicator.h"

int main(int argc, char** argv)
{
	MPIEnvironment mpi(argc, argv);
	const Communicator& world = mpi.World();
	CurrentCommunicator current(world);
	::testing::InitGoogleTest(&argc, argv);
	if (!world.IsRoot())
	{
		::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
		delete listeners.Release(listeners.default_result_printer());
		cout.rdbuf(nullptr);
	}
	int result = RUN_ALL_TESTS();
	return world.Max(result);
}
