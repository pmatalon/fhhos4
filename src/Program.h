#pragma once
#include "Program/Program_Diffusion_DG.h"
#include "Program/Program_Diffusion_HHO.h"
#include "Program/Program_BiHarmonic_HHO.h"
using namespace std;

class Program
{
public:
	Program() {}
	virtual void Start(ProgramArguments& args) = 0;
	virtual ~Program() {}
};

template <int Dim>
class ProgramDim : public Program
{
public:
	ProgramDim() : Program() {}

	void Start(ProgramArguments& args)
	{
		Utils::ProgramArgs = args;

		Timer totalTimer;
		totalTimer.Start();

#ifdef SMALL_INDEX
		cout << "Index type: int" << endl;
#else
		cout << "Index type: size_t" << endl;
#endif
		cout << "Shared memory parallelism: " << (BaseParallelLoop::GetDefaultNThreads() == 1 ? "sequential execution" : to_string(BaseParallelLoop::GetDefaultNThreads()) + " threads") << endl;

		Mesh<Dim>::SetDirectories();
		GMSHMesh<Dim>::GMSHLogEnabled = args.Actions.GMSHLogEnabled;
		GMSHMesh<Dim>::UseCache = args.Actions.UseCache;

		if (args.Problem.Equation == EquationType::Diffusion)
		{
			if (args.Discretization.Method.compare("dg") == 0)
				Program_Diffusion_DG<Dim>::Execute(args);
			else if (args.Discretization.Method.compare("hho") == 0)
				Program_Diffusion_HHO<Dim>::Execute(args);
			else
				Utils::FatalError("Unknown or unmanaged discretization for diffusion problem. Check arguments -pb and -discr.");
		}
		else if (args.Problem.Equation == EquationType::BiHarmonic)
		{
			if (args.Discretization.Method.compare("hho") == 0)
				Program_BiHarmonic_HHO<Dim>::Execute(args);
			else
				Utils::FatalError("Unknown or unmanaged discretization for bi-harmonic problem. Check arguments -pb and -discr.");
		}
		else
			Utils::FatalError("Unknown problem. Check argument -pb.");

		totalTimer.Stop();
		cout << endl << "Total time: CPU = " << totalTimer.CPU() << ", elapsed = " << totalTimer.Elapsed() << endl;
	}

};
