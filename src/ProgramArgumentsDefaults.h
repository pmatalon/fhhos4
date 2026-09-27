#pragma once
#include "Program.h"
using namespace std;

// Prints the message as a fatal argument error and terminates the process.
// Shared by main()'s CLI parsing and by ApplyProgramArgumentDefaults() below.
inline void argument_error(string msg)
{
	cout << Utils::BeginRed << "Argument error: " << msg << Utils::EndColor << endl;
	cout << "------------------------- FAILURE -------------------------" << endl;
	exit(EXIT_FAILURE);
}

// Resolves every sentinel/"default" value left in a ProgramArguments (Dimension == -1,
// SolverCode/Mesher/MeshCode == "default", etc.) into a concrete, consistent configuration,
// exactly as main()'s CLI parsing does after getopt. Callers that build a ProgramArguments
// directly (e.g. tests) can call this instead of going through the CLI to get the same
// defaulting behavior.
//
// The boolean flags mirror whether the corresponding CLI option was explicitly passed by the
// user (main() tracks them while parsing getopt); left at their default of `true` (= "not
// explicitly set"), the same default-inference cascade as the CLI's own defaults applies.
inline void ApplyProgramArgumentDefaults(ProgramArguments& args,
	bool defaultRelativeCellPolyDegree = true,
	bool defaultTol2 = true,
	bool defaultCycle = true,
	bool defaultCoarseOperator = true,
	bool defaultCoarseSolver = true)
{
	//------------------------------------------//
	//                 Problem                  //
	//------------------------------------------//

	if (args.Problem.Equation == EquationType::BiHarmonic && defaultRelativeCellPolyDegree)
		args.Discretization.RelativeCellPolyDegree = 1;

	// Dimension
	if (args.Problem.Dimension == -1)
	{
		if (args.Problem.GeoCode.compare("segment") == 0)
			args.Problem.Dimension = 1;
		else if (args.Problem.GeoCode.compare("square") == 0 || args.Problem.GeoCode.compare("square4quadrants") == 0 || args.Problem.GeoCode.compare("L_shape") == 0)
			args.Problem.Dimension = 2;
		else if (args.Problem.GeoCode.compare("cube") == 0)
			args.Problem.Dimension = 3;
#ifdef GMSH_ENABLED
		else
		{
#ifdef ENABLE_2D
			Mesh<2>::SetDirectories();
			args.Problem.Dimension = GMSHMesh<2>::GetDimension(args.Problem.GeoCode);
#elif defined ENABLE_3D
			Mesh<3>::SetDirectories();
			args.Problem.Dimension = GMSHMesh<3>::GetDimension(args.Problem.GeoCode);
#endif // ENABLE_2D
		}
#else
		else
			argument_error("Unknown geometry.");
#endif // GMSH_ENABLED
	}

	// Test case
	if (args.Problem.TestCaseCode.compare("") == 0)
	{
		if (Utils::IsPredefinedGeometry(args.Problem.GeoCode))
			args.Problem.TestCaseCode = args.Problem.GeoCode;
		else
			args.Problem.TestCaseCode = FileSystem::FileNameWithoutExtension(args.Problem.GeoCode);
	}

	// Heterogeneity
	if (args.Problem.GeoCode.compare("square") == 0 && args.Problem.HeterogeneityRatio != 1)
		Utils::Warning("The geometry 'square' has only one physical part: -heterog argument is ignored. Use 'square4quadrants' instead to run an heterogeneous problem.");

	// Boundary conditions
	if (args.Problem.BCCode.compare("") == 0)
	{
		if (args.Problem.TestCaseCode.compare("fullneumann") == 0)
			args.Problem.BCCode = "n2";
		else
			args.Problem.BCCode = "d";
	}

	//------------------------------------------//
	//                  Mesher                  //
	//------------------------------------------//

	if (args.Discretization.Mesher.compare("default") == 0)
	{
		if (args.Problem.GeoCode.compare("cube") == 0 && args.Discretization.MeshCode.compare("cart") == 0)
			args.Discretization.Mesher = "inhouse";
		else
#ifdef GMSH_ENABLED
			args.Discretization.Mesher = "gmsh";
#else
			args.Discretization.Mesher = "inhouse";
#endif
	}


#ifndef GMSH_ENABLED
	if (args.Discretization.Mesher.compare("gmsh") == 0)
		argument_error("GMSH is disabled. Recompile with the cmake option -DENABLE_GMSH=ON to use GMSH meshes, or choose another argument for -mesh.");
#endif // GMSH_ENABLED

	//------------------------------------------//
	//                   Mesh                   //
	//------------------------------------------//

	if (args.Discretization.MeshCode.compare("default") == 0)
	{
		if (args.Problem.Dimension == 1)
			args.Discretization.MeshCode = "cart";
		else if (args.Problem.Dimension == 2)
			args.Discretization.MeshCode = args.Discretization.Mesher.compare("inhouse") == 0 ? "stri" : "tri";
		else if (args.Problem.Dimension == 3)
			args.Discretization.MeshCode = args.Discretization.Mesher.compare("inhouse") == 0 ? "stetra" : "tetra";
	}

	if ((args.Discretization.MeshCode.compare("tri") == 0 || args.Discretization.MeshCode.compare("stri") == 0) && args.Problem.Dimension != 2)
		argument_error("Triangular mesh in only available in 2D.");

	if (args.Discretization.MeshCode.compare("quad") == 0 && args.Problem.Dimension != 2)
		argument_error("Quadrilateral mesh in only available in 2D.");

	if ((args.Discretization.MeshCode.compare("tetra") == 0 || args.Discretization.MeshCode.compare("stetra") == 0) && args.Problem.Dimension != 3)
		argument_error("Tetrahedral mesh in only available in 3D.");

#ifndef CGAL_ENABLED
	if (args.Discretization.MeshCode.compare("poly") == 0)
		Utils::FatalError("CGAL must be enabled to use polygonal meshes. Recompile the program with cmake option -DENABLE_CGAL=On.");
#endif

	//------------------------------------------//
	//              Discretization              //
	//------------------------------------------//

	if (args.Problem.Dimension > 1 && args.Discretization.Method.compare("dg") == 0 && args.Discretization.PolyDegree == 0)
		argument_error("In 2D/3D, DG is not a convergent scheme for p = 0.");

	if (args.Discretization.Method.compare("dg") == 0 && args.Problem.BCCode.compare("d") != 0)
		argument_error("In DG, only Dirichlet conditions are implemented.");

	if (args.Discretization.Method.compare("dg") == 0 && args.Problem.AnisotropyRatio != 1)
		argument_error("In DG, anisotropy is not implemented.");

	if (args.Problem.Dimension == 1 && args.Discretization.Method.compare("hho") == 0 && args.Discretization.PolyDegree != 1)
		argument_error("HHO in 1D only exists for k = 0.");

	if (args.Discretization.Method.compare("hho") == 0 && args.Discretization.PolyDegree == 0)
		argument_error("HHO does not exist with p = 0. Linear approximation at least (p >= 1).");

	if (args.Discretization.Method.compare("hho") == 0 && args.Discretization.Stabilization.compare("hdg") == 0 && args.Discretization.RelativeCellPolyDegree < 1)
	{
		if (args.Discretization.PolyDegree == 1)
			Utils::Warning("HHO(k=0) is not convergent with the 'hdg' stabilization. With '-stab hdg', you should use '-kc 1'.");
		else
			Utils::Warning("HHO loses one degree of approximation with the 'hdg' stabilization. With '-stab hdg', you should use '-kc 1'.");
	}

	// Elem polynomial bases
	if (args.Discretization.Method.compare("fem") == 0)
	{
		if (!args.Discretization.ElemBasisCode.empty() && args.Discretization.ElemBasisCode.compare("lagrange") != 0)
			argument_error("In FEM, only Lagrange basis is implemented.");
		if (args.Discretization.PolyDegree != 1)
			argument_error("In FEM, only p=1 is implemented.");
		args.Discretization.ElemBasisCode = "lagrange";
		args.Discretization.OrthogonalizeElemBasesCode = 0;
		args.Discretization.PolyDegree = 1;
	}
	else if (args.Discretization.ElemBasisCode.empty())
	{
		if (args.Discretization.OrthogonalizeElemBasesCode == -1)
			args.Discretization.OrthogonalizeElemBasesCode = args.Discretization.MeshCode.compare("cart") == 0 ? 0 : 1;
		args.Discretization.ElemBasisCode = args.Discretization.OrthogonalizeElemBasesCode > 0 ? "monomials" : "legendre";
	}
	else if (args.Discretization.OrthogonalizeElemBasesCode == -1)
	{
		if (args.Discretization.ElemBasisCode.compare("legendre") == 0 && args.Discretization.MeshCode.compare("cart") == 0)
			args.Discretization.OrthogonalizeElemBasesCode = 0;
		else
			args.Discretization.OrthogonalizeElemBasesCode = 1;
	}

	// Face polynomial bases
	if (args.Discretization.FaceBasisCode.empty())
	{
		if (args.Discretization.OrthogonalizeFaceBasesCode == -1)
		{
			if (args.Problem.Dimension <= 2 || args.Discretization.MeshCode.compare("cart") == 0)
				args.Discretization.OrthogonalizeFaceBasesCode = 0;
			else
				args.Discretization.OrthogonalizeFaceBasesCode = 1;
		}
		args.Discretization.FaceBasisCode = args.Discretization.OrthogonalizeFaceBasesCode > 0 ? "monomials" : "legendre";
	}
	else if (args.Discretization.OrthogonalizeFaceBasesCode == -1)
	{
		if (args.Discretization.FaceBasisCode.compare("legendre") == 0 && args.Discretization.MeshCode.compare("cart") == 0)
			args.Discretization.OrthogonalizeFaceBasesCode = 0;
		else
			args.Discretization.OrthogonalizeFaceBasesCode = 1;
	}

	//------------------------------------------//
	//                  Solver                  //
	//------------------------------------------//

	if (args.Solver.SolverCode.compare("default") == 0)
	{
		if (args.Discretization.Method.compare("hho") == 0 && args.Discretization.StaticCondensation && args.Problem.Dimension > 1)
		{
			args.Solver.SolverCode = "mg";
			if ((args.Solver.MG.GMG_H_Prolong == GMG_H_Prolongation::Wildey || args.Solver.MG.GMG_H_Prolong == GMG_H_Prolongation::FaceInject) && !args.Solver.MG.UseGalerkinOperator)
			{
				Utils::Warning("The multigrid with prolongation code " + to_string((unsigned)args.Solver.MG.GMG_H_Prolong) + " requires the Galerkin operator. Option -g 0 ignored.");
				args.Solver.MG.UseGalerkinOperator = true;
			}
		}
		else if ((args.Problem.Dimension == 2 && args.Discretization.N < 64) || (args.Problem.Dimension == 3 && args.Discretization.N < 16))
			args.Solver.SolverCode = "lu";
		else
			args.Solver.SolverCode = "eigencg";
	}
	// Retrocompatibility
	else if (Utils::StartsWith(args.Solver.SolverCode, "cg") && args.Solver.SolverCode.length() > 2)
	{
		args.Solver.PreconditionerCode = args.Solver.SolverCode.substr(2, args.Solver.SolverCode.length() - 2);
		args.Solver.SolverCode = "cg";
	}
	else if (Utils::StartsWith(args.Solver.SolverCode, "fcg") && args.Solver.SolverCode.length() > 3)
	{
		args.Solver.PreconditionerCode = args.Solver.SolverCode.substr(3, args.Solver.SolverCode.length() - 3);
		args.Solver.SolverCode = "fcg";
	}

#ifndef AGMG_ENABLED
	if (args.Solver.SolverCode.compare("agmg") == 0)
		argument_error("AGMG is disabled. Recompile with the cmake option -DENABLE_AGMG=ON, or choose another solver.");
#endif // AGMG_ENABLED

	if (defaultTol2)
		args.Solver.Tolerance2 = args.Solver.Tolerance;

	//------------------------------------------//
	//                Multigrid                 //
	//------------------------------------------//

	if (args.Solver.SolverCode.compare("mg") == 0 || args.Solver.PreconditionerCode.compare("mg") == 0 || args.Solver.SolverCode.compare("p_mg") == 0 || args.Solver.PreconditionerCode.compare("p_mg") == 0)
	{
		if (args.Discretization.Method.compare("dg") == 0)
			argument_error("Multigrid only applicable on HHO discretization.");

		if (!args.Discretization.StaticCondensation)
			argument_error("Multigrid only applicable if the static condensation is enabled.");

		args.Solver.MG.GMG_H_Prolong = static_cast<GMG_H_Prolongation>(args.Solver.MG.ProlongationCode);

		if (args.Solver.MG.GMG_H_Prolong == GMG_H_Prolongation::Wildey && !args.Solver.MG.UseGalerkinOperator)
			argument_error("To use the prolongationCode " + to_string((unsigned)GMG_H_Prolongation::Wildey) + ", you must also use the Galerkin operator. To do so, add option -g 1.");

		if (args.Solver.MG.H_CS == H_CoarsStgy::FaceCoarsening && !args.Solver.MG.UseGalerkinOperator)
			argument_error("To use the face coarsening, you must also use the Galerkin operator. To do so, add option -g 1.");

		if (args.Solver.MG.H_CS == H_CoarsStgy::IndependentRemeshing &&
			Utils::RequiresNestedHierarchy(args.Solver.MG.GMG_H_Prolong) &&
			args.Solver.MG.GMG_H_Prolong != GMG_H_Prolongation::Default)
			argument_error("The coarsening by independent remeshing is only applicable with the non-nested versions of the multigrid (-prolong " + to_string((unsigned)GMG_H_Prolongation::CellInterp_ExactL2proj_Trace) + ", " + to_string((unsigned)GMG_H_Prolongation::CellInterp_ApproxL2proj_Trace) + " or " + to_string((unsigned)GMG_H_Prolongation::CellInterp_FinerApproxL2proj_Trace) + ").");

		if (args.Solver.SolverCode.compare("p_mg") == 0 && defaultCoarseSolver)
			args.Solver.MG.CoarseSolverCode = "mg";

		if (args.Solver.MG.H_CS == H_CoarsStgy::None)
		{
			if (args.Problem.Dimension < 3)
			{
				if (args.Solver.MG.GMG_H_Prolong == GMG_H_Prolongation::Default)
				{
					if (args.Discretization.Mesher.compare("inhouse") == 0)
						args.Solver.MG.H_CS = H_CoarsStgy::StandardCoarsening;
					else
						args.Solver.MG.H_CS = H_CoarsStgy::IndependentRemeshing;
				}
				else if (Utils::RequiresNestedHierarchy(args.Solver.MG.GMG_H_Prolong))
					args.Solver.MG.H_CS = H_CoarsStgy::GMSHSplittingRefinement;
				else
					args.Solver.MG.H_CS = H_CoarsStgy::IndependentRemeshing;
			}
			else
			{
				if (args.Discretization.Mesher.compare("inhouse") == 0 && args.Discretization.MeshCode.compare("tetra") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::BeyRefinement;
				else if (args.Discretization.Mesher.compare("inhouse") == 0 && args.Discretization.MeshCode.compare("cart") == 0)
					args.Solver.MG.H_CS = H_CoarsStgy::StandardCoarsening;
				else
					args.Solver.MG.H_CS = H_CoarsStgy::IndependentRemeshing;
			}
		}

		if (args.Solver.MG.GMG_H_Prolong == GMG_H_Prolongation::Default)
		{
			if (Utils::BuildsNestedMeshHierarchy(args.Solver.MG.H_CS))
				args.Solver.MG.GMG_H_Prolong = GMG_H_Prolongation::CellInterp_Trace;
			else
				args.Solver.MG.GMG_H_Prolong = GMG_H_Prolongation::CellInterp_FinerApproxL2proj_Trace;
		}

		if (defaultCycle)
		{
			args.Solver.MG.PreSmoothingIterations = 0;
			args.Solver.MG.PostSmoothingIterations = args.Problem.Dimension < 3 ? 3 : 6;
		}
	}

	if (args.Solver.SolverCode.compare("uamg") == 0 || args.Solver.PreconditionerCode.compare("uamg") == 0)
	{
		if (args.Discretization.Method.compare("dg") == 0)
			argument_error("Multigrid only applicable on HHO discretization.");

		if (args.Solver.MG.ProlongationCode == 0)
			args.Solver.MG.UAMGMultigridProlong = UAMGProlongation::ChainedCoarseningProlongations;
		else
			args.Solver.MG.UAMGMultigridProlong = static_cast<UAMGProlongation>(args.Solver.MG.ProlongationCode);
		args.Solver.MG.UAMGFaceProlong = args.Solver.MG.FaceProlongationCode == 0 ? UAMGFaceProlongation::BoundaryAggregatesInteriorAverage : static_cast<UAMGFaceProlongation>(args.Solver.MG.FaceProlongationCode);
		args.Solver.MG.UAMGCoarseningProlong = args.Solver.MG.CoarseningProlongationCode == 0 ? UAMGProlongation::ReconstructSmoothedTraceOrInject : static_cast<UAMGProlongation>(args.Solver.MG.CoarseningProlongationCode);

		if (args.Solver.MG.H_CS == H_CoarsStgy::None)
			args.Solver.MG.H_CS = H_CoarsStgy::MultiplePairwiseAggregation;

		if ((args.Solver.MG.H_CS == H_CoarsStgy::MultiplePairwiseAggregation ||
			 args.Solver.MG.H_CS == H_CoarsStgy::MultipleAgglomerationCoarseningByFaceNeighbours)
			&& args.Solver.MG.CoarseningFactor == 0)
			args.Solver.MG.CoarseningFactor = 3.8;

		if (defaultCoarseOperator)
			args.Solver.MG.UseGalerkinOperator = true;

		if (defaultCycle)
			args.Solver.MG.CycleLetter = 'K';
	}

	if (args.Solver.SolverCode.compare("aggregamg") == 0 || args.Solver.PreconditionerCode.compare("aggregamg") == 0)
	{
		if (!defaultCoarseOperator && !args.Solver.MG.UseGalerkinOperator)
			Utils::Warning("AggregAMG uses the Galerkin operator. Argument -g 0 ignored.");

		if (args.Solver.MG.H_CS == H_CoarsStgy::None)
			args.Solver.MG.H_CS = H_CoarsStgy::DoublePairwiseAggregation;

		if (defaultCycle)
			args.Solver.MG.CycleLetter = 'K';
	}


	if (args.Solver.MG.H_CS == H_CoarsStgy::MultiplePairwiseAggregation && args.Solver.MG.CoarseningFactor == 0)
		args.Solver.MG.CoarseningFactor = 3.5;
	else if (args.Solver.MG.H_CS == H_CoarsStgy::IndependentRemeshing && args.Solver.MG.CoarseningFactor == 0)
		args.Solver.MG.CoarseningFactor = 2;
}
