#include "Program.h"
#include "Program/Program_Diffusion_DG_Impl.h"
#include "Program/Program_Diffusion_HHO_Impl.h"
#include "Program/Program_Diffusion_FEM_Impl.h"
#include "Program/Program_BiHarmonic_HHO_Impl.h"
#include "Program/Program_BiHarmonic_FEM_Impl.h"
#include "Program/Program_BiHarmonicDD_HHO_Impl.h"
#include "Mesher/GMSH/GMSHMesh.h"
#include "Utils/Timer.h"
#ifdef CGAL_ENABLED
#include "Geometry/CGALWrapper.h"
#endif

template <int Dim>
void ProgramDim<Dim>::InitGlobalState(const ProgramArguments& args)
{
	Utils::ProgramArgs = args;
	Mesh<Dim>::SetDirectories();
	GMSHMesh<Dim>::GMSHLogEnabled = args.Actions.GMSHLogEnabled;
	GMSHMesh<Dim>::UseCache = args.Actions.UseCache;
}

template <int Dim>
void ProgramDim<Dim>::Start(ProgramArguments& args)
{
#ifdef CGAL_ENABLED
	CGALWrapper::Configure();
#endif
	InitGlobalState(args);

	Timer totalTimer;
	totalTimer.Start();

#ifdef SMALL_INDEX
	cout << "Index type: int" << endl;
#else
	cout << "Index type: size_t" << endl;
#endif
	cout << "Shared memory parallelism: " << (BaseParallelLoop::GetDefaultNThreads() == 1 ? "sequential execution" : to_string(BaseParallelLoop::GetDefaultNThreads()) + " threads") << endl;
	cout << endl;

	if (args.Problem.Equation == EquationType::Diffusion)
	{
		if (args.Discretization.Method.compare("dg") == 0)
			Program_Diffusion_DG<Dim>::Execute(args);
		else if (args.Discretization.Method.compare("hho") == 0)
			Program_Diffusion_HHO<Dim>::Execute(args); 
		else if (args.Discretization.Method.compare("fem") == 0)
			Program_Diffusion_FEM<Dim>::Execute(args);
		else
			Utils::FatalError("Unknown or unmanaged discretization for diffusion problem. Check arguments -pb and -discr.");
	}
	else if (args.Problem.Equation == EquationType::BiHarmonic)
	{
		if (args.Discretization.Method.compare("hho") == 0)
			Program_BiHarmonic_HHO<Dim>::Execute(args);
		else if (args.Discretization.Method.compare("fem") == 0)
			Program_BiHarmonic_FEM<Dim>::Execute(args);
		else
			Utils::FatalError("Unknown or unmanaged discretization for bi-harmonic problem. Check arguments -pb and -discr.");
	}
	else if (args.Problem.Equation == EquationType::BiHarmonicDD)
	{
		if (args.Discretization.Method.compare("hho") == 0)
			Program_BiHarmonicDD_HHO<Dim>::Execute(args);
		else
			Utils::FatalError("Unknown or unmanaged discretization for bi-harmonic problem. Check arguments -pb and -discr.");
	}
	else
		Utils::FatalError("Unknown problem. Check argument -pb.");

	totalTimer.Stop();
	cout << endl << "Total time: CPU = " << totalTimer.CPU() << ", elapsed = " << totalTimer.Elapsed() << endl;
}

// The programs are only compiled here: the shared infrastructure they instantiate (meshes,
// solvers, bases...) dominates the compilation time, and would be compiled again in every
// translation unit instantiating a program.
#ifdef ENABLE_1D
template class ProgramDim<1>;
template class Program_BiHarmonicDD_HHO<1>;
template class Program_BiHarmonic_FEM<1>;
template class Program_BiHarmonic_HHO<1>;
template class Program_Diffusion_DG<1>;
template class Program_Diffusion_FEM<1>;
template class Program_Diffusion_HHO<1>;
#endif // ENABLE_1D

#ifdef ENABLE_2D
template class ProgramDim<2>;
template class Program_BiHarmonicDD_HHO<2>;
template class Program_BiHarmonic_FEM<2>;
template class Program_BiHarmonic_HHO<2>;
template class Program_Diffusion_DG<2>;
template class Program_Diffusion_FEM<2>;
template class Program_Diffusion_HHO<2>;
#endif // ENABLE_2D

#ifdef ENABLE_3D
template class ProgramDim<3>;
template class Program_BiHarmonicDD_HHO<3>;
template class Program_BiHarmonic_FEM<3>;
template class Program_BiHarmonic_HHO<3>;
template class Program_Diffusion_DG<3>;
template class Program_Diffusion_FEM<3>;
template class Program_Diffusion_HHO<3>;
#endif // ENABLE_3D
