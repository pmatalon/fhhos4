#pragma once
#include "../ProgramArguments.h"
#include "../TestCases/TestCase.h"

#ifdef ENABLE_1D
	#include "InHouse/UniformMesh1D.h"
#endif // ENABLE_1D

#ifdef ENABLE_2D
	#include "InHouse/Square_CartesianMesh.h"
	#include "InHouse/Square_CartesianPolygonalMesh.h"
	#include "InHouse/Square_TriangularMesh.h"
	#include "InHouse/Square_QuadrilateralMesh.h"
	#ifdef CGAL_ENABLED
		#include "InHouse/Square_QuadrilateralAsPolygonalMesh.h"
	#endif // CGAL_ENABLED
	#ifdef GMSH_ENABLED
		#include "GMSH/Square_GMSHCartesianMesh.h"
		#include "GMSH/Square_GMSHTriangularMesh.h"
		#include "GMSH/Square_GMSHUnstructTriangularMesh.h"
		#include "GMSH/Square_GMSHQuadrilateralMesh.h"

		#include "GMSH/Square4quadrants_GMSHCartesianMesh.h"
		#include "GMSH/Square4quadrants_GMSHTriangularMesh.h"
		#include "GMSH/Square4quadrants_GMSHUnstructTriangularMesh.h"
		#include "GMSH/Square4quadrants_GMSHQuadrilateralMesh.h"
	#endif // GMSH_ENABLED
#endif // ENABLE_2D

#ifdef ENABLE_3D
	#include "InHouse/Cube_CartesianMesh.h"
	#include "InHouse/Cube_CartesianTetrahedralMesh.h"
	#ifdef GMSH_ENABLED
		#include "GMSH/Cube_GMSHTetrahedralMesh.h"
		#include "GMSH/Cube_GMSHCartesianMesh.h"
	#endif // GMSH_ENABLED
#endif // ENABLE_3D

using namespace std;

template <int Dim>
class MeshFactory
{
public:
	static Mesh<Dim>* BuildMesh(ProgramArguments& args, TestCase<Dim>* testCase) { return nullptr; }
private:
	static PolyhedralMesh<Dim>* BuildPolyhedralMesh(PolyhedralMesh<Dim>* meshToAggregate, FaceCoarseningStrategy faceCoarseningStgy, FaceCollapsing bdryFaceCollapsing) { return nullptr; }
};

#ifdef ENABLE_1D
template <>
Mesh<1>* MeshFactory<1>::BuildMesh(ProgramArguments& args, TestCase<1>* testCase)
{
	return new UniformMesh1D(args.Discretization.N);
}
#endif // ENABLE_1D

#ifdef ENABLE_2D

template <>
PolyhedralMesh<2>* MeshFactory<2>::BuildPolyhedralMesh(PolyhedralMesh<2>* mesh, FaceCoarseningStrategy faceCoarseningStgy, FaceCollapsing bdryFaceCollapsing)
{
	cout << "Building polygonal mesh by agglomeration:" << endl;
	cout << "\t" << "Coarsening strategy     : agglomeration by face neighbours" << endl;

	cout << "\t" << "Face coarsening strategy: ";
	if (faceCoarseningStgy == FaceCoarseningStrategy::None)
		cout << "no coarsening" << endl;
	else if (faceCoarseningStgy == FaceCoarseningStrategy::InterfaceCollapsing)
		cout << "interface collapsing [-fcs c]" << endl;
	else
		cout << "unknown" << endl;

	cout << "\t" << "Boundary face collapsing: ";
	if (bdryFaceCollapsing == FaceCollapsing::Disabled)
		cout << "disabled [-bfc d]" << endl;
	else if (bdryFaceCollapsing == FaceCollapsing::OnlyCollinear)
		cout << "collinear only [-bfc c]" << endl;
	else if (bdryFaceCollapsing == FaceCollapsing::ByPairs)
		cout << "by pairs [-bfc p]" << endl;
	else if (bdryFaceCollapsing == FaceCollapsing::Max)
		cout << "maximum [-bfc m]" << endl;
	else
		cout << "unknown" << endl;

	mesh->CoarsenMesh(H_CoarsStgy::AgglomerationCoarseningByFaceNeighbours, faceCoarseningStgy, bdryFaceCollapsing, 0);
	PolyhedralMesh<2>* polyMesh = static_cast<PolyhedralMesh<2>*>(mesh->CoarseMesh);
	mesh->CoarseMesh = nullptr;
	polyMesh->FineMesh = nullptr;
	polyMesh->Vertices = std::move(mesh->Vertices);
	polyMesh->ClearMeshVertexConnections();
	polyMesh->UpdateMeshVertexConnections();
	return polyMesh;
}



template <>
Mesh<2>* MeshFactory<2>::BuildMesh(ProgramArguments& args, TestCase<2>* testCase)
{
	string geoCode = args.Problem.GeoCode;
	string mesher = args.Discretization.Mesher;
	BigNumber n = args.Discretization.N;
	BigNumber nx = args.Discretization.N;
	BigNumber ny = args.Discretization.Ny == 0 ? args.Discretization.N : args.Discretization.Ny;
	string meshCode = args.Discretization.MeshCode;
	double stretch = args.Discretization.Stretch;

	H_CoarsStgy refinementStgy = Utils::IsRefinementStrategy(args.Solver.MG.H_CS) ? args.Solver.MG.H_CS : H_CoarsStgy::GMSHSplittingRefinement;

	if (refinementStgy == H_CoarsStgy::BeyRefinement)
		Utils::FatalError("Bey's refinement method is only applicable to 3D tetrahedral meshes.");

	Mesh<2>* fineMesh = nullptr;
	//-------------------//
	//       Square      //
	//-------------------//
	if (geoCode.compare("square") == 0)
	{
		if (mesher.compare("inhouse") == 0)
		{
			if (Utils::IsRefinementStrategy(args.Solver.MG.H_CS))
				Utils::FatalError("Unmanaged refinement strategy.");

			if (meshCode.compare("cart") == 0)
				fineMesh = new Square_CartesianMesh(nx, ny);
			else if (meshCode.compare("cart-poly") == 0)
				fineMesh = new Square_CartesianPolygonalMesh(nx, ny);
			else if (meshCode.compare("stri") == 0)
				fineMesh = new Square_TriangularMesh(nx, ny);
			else if (meshCode.compare("quad") == 0)
				fineMesh = new Square_QuadrilateralMesh(nx, ny, stretch);
#ifdef CGAL_ENABLED
			else if (meshCode.compare("quad-poly") == 0)
				fineMesh = new Square_QuadrilateralAsPolygonalMesh(nx, ny, stretch);
#endif
			else if (meshCode.compare("tri") == 0)
				Utils::FatalError("The in-house mesher does not build unstructured meshes. Use '-mesher gmsh' or '-mesh stri' instead.");
			else
				Utils::FatalError("The requested mesh is not managed with this geometry.");
		}
#ifdef GMSH_ENABLED
		else if (mesher.compare("gmsh") == 0)
		{
			if (meshCode.compare("cart") == 0)
			{
				if (args.Solver.MG.H_CS == H_CoarsStgy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square_GMSHCartesianMesh(n <= 16 ? 2 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(nx * ny, refinementStgy);
				}
				else
					fineMesh = new Square_GMSHCartesianMesh(n);
			}
			else if (meshCode.compare("stri") == 0)
			{
				if (args.Solver.MG.H_CS == H_CoarsStgy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square_GMSHTriangularMesh(n <= 16 ? 2 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(2 * nx * ny, refinementStgy);
				}
				else
					fineMesh = new Square_GMSHTriangularMesh(n);
			}
			else if (meshCode.compare("tri") == 0)
			{
				if (args.Solver.MG.H_CS == H_CoarsStgy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square_GMSHUnstructTriangularMesh(n <= 16 ? 2 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(2 * nx * ny, refinementStgy);
				}
				else
					fineMesh = new Square_GMSHUnstructTriangularMesh(n);
			}
			else if (meshCode.compare("quad") == 0)
			{
				if (args.Solver.MG.H_CS == H_CoarsStgy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square_GMSHQuadrilateralMesh(n <= 16 ? 2 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(nx * ny, refinementStgy);
				}
				else
					fineMesh = new Square_GMSHQuadrilateralMesh(n);
			}
			else if (meshCode.compare("poly") == 0)
			{
				Square_GMSHUnstructTriangularMesh* triMesh = new Square_GMSHUnstructTriangularMesh(n);
				fineMesh = BuildPolyhedralMesh(triMesh, args.Discretization.PolyMeshFaceCoarseningStgy, args.Discretization.PolyMeshBoundaryFaceCollapsing);
			}
			else
				Utils::FatalError("The requested mesh is not managed with this geometry.");
		}
#endif // GMSH_ENABLED
		else
			Utils::FatalError("Unknown mesher. Check -mesher argument.");
	}
	//-------------------------------//
	//       Square 4 quadrants      //
	//-------------------------------//
	else if (geoCode.compare("square4quadrants") == 0)
	{
		bool with4quadrants = true;
		if (mesher.compare("inhouse") == 0)
		{
			if (Utils::IsRefinementStrategy(args.Solver.MG.H_CS))
				Utils::FatalError("Unmanaged refinement strategy.");

			if (meshCode.compare("cart") == 0)
				fineMesh = new Square_CartesianMesh(nx, ny, with4quadrants);
			else if (meshCode.compare("cart-poly") == 0)
				fineMesh = new Square_CartesianPolygonalMesh(nx, ny, with4quadrants);
			else if (meshCode.compare("stri") == 0)
				fineMesh = new Square_TriangularMesh(nx, ny, with4quadrants);
			else if (meshCode.compare("tri") == 0)
				Utils::FatalError("The in-house mesher does not build unstructured meshes. Use '-mesh stri' or '-mesher gmsh' instead.");
			else
				Utils::FatalError("The requested mesh is not managed with this geometry.");
		}
#ifdef GMSH_ENABLED
		else if (mesher.compare("gmsh") == 0)
		{
			if (meshCode.compare("cart") == 0)
			{
				if (args.Solver.MG.H_CS == H_CoarsStgy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square4quadrants_GMSHCartesianMesh(n <= 16 ? 4 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(nx * ny, refinementStgy);
				}
				else
					fineMesh = new Square4quadrants_GMSHCartesianMesh(n);
			}
			else if (meshCode.compare("tri") == 0)
			{
				if (args.Solver.MG.H_CS == H_CoarsStgy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square4quadrants_GMSHUnstructTriangularMesh(n <= 16 ? 4 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(2 * nx * ny, refinementStgy);
				}
				else
					fineMesh = new Square4quadrants_GMSHUnstructTriangularMesh(n);
			}
			else if (meshCode.compare("stri") == 0)
			{
				if (args.Solver.MG.H_CS == H_CoarsStgy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square4quadrants_GMSHTriangularMesh(n <= 16 ? 4 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(2 * nx * ny, refinementStgy);
				}
				else
					fineMesh = new Square4quadrants_GMSHTriangularMesh(n);
			}
			else if (meshCode.compare("quad") == 0)
			{
				if (args.Solver.MG.H_CS == H_CoarsStgy::GMSHSplittingRefinement)
				{
					Mesh<2>* coarseMesh = new Square4quadrants_GMSHQuadrilateralMesh(n <= 16 ? 4 : 16);
					fineMesh = coarseMesh->RefineUntilNElements(nx * ny, refinementStgy);
				}
				else
					fineMesh = new Square4quadrants_GMSHQuadrilateralMesh(n);
			}
			else
				Utils::FatalError("The requested mesh is not managed with this geometry.");
		}
#endif // GMSH_ENABLED
		else
			Utils::FatalError("Unknown mesher.");
	}
	//----------------------//
	//       GMSH file      //
	//----------------------//
	else
	{
		if (mesher.compare("inhouse") == 0)
			Utils::FatalError("The geometry is imported from a GMSH file, the mesher should be 'gmsh'. Use '-mesher gmsh' instead.");

#ifdef GMSH_ENABLED
		string filePath = geoCode;

		if (meshCode.compare("poly") == 0)
		{
			GMSHMesh<2>* mesh = new GMSHMesh<2>(testCase, filePath, n);
			fineMesh = BuildPolyhedralMesh(mesh, args.Discretization.PolyMeshFaceCoarseningStgy, args.Discretization.PolyMeshBoundaryFaceCollapsing);
		}
		else if (args.Solver.MG.H_CS == H_CoarsStgy::GMSHSplittingRefinement)
		{
			Mesh<2>* coarseMesh = new GMSHMesh<2>(testCase, filePath, args.Solver.MG.CoarseN);
			fineMesh = coarseMesh->RefineUntilNElements(2 * nx * ny, refinementStgy);
		}
		else
			fineMesh = new GMSHMesh<2>(testCase, filePath, n);
#endif // GMSH_ENABLED
	}

	if (!Utils::IsRefinementStrategy(args.Solver.MG.H_CS) && fineMesh->CoarseMesh)
		fineMesh->DeleteCoarseMeshes();

	return fineMesh;
}
#endif // ENABLE_2D



#ifdef ENABLE_3D

template <>
Mesh<3>* MeshFactory<3>::BuildMesh(ProgramArguments& args, TestCase<3>* testCase)
{
	string geoCode = args.Problem.GeoCode;
	string mesher = args.Discretization.Mesher;
	BigNumber n = args.Discretization.N;
	BigNumber nx = args.Discretization.N;
	BigNumber ny = args.Discretization.Ny == 0 ? args.Discretization.N : args.Discretization.Ny;
	BigNumber nz = args.Discretization.Nz == 0 ? args.Discretization.N : args.Discretization.Nz;
	string meshCode = args.Discretization.MeshCode;
	H_CoarsStgy refinementStgy = args.Solver.MG.H_CS;

	Mesh<3>* fineMesh = nullptr;
	//------------------------//
	//          Cube          //
	//------------------------//
	if (geoCode.compare("cube") == 0)
	{
		if (mesher.compare("inhouse") == 0)
		{
			if (meshCode.compare("cart") == 0)
			{
				fineMesh = new Cube_CartesianMesh(nx, ny, nz);

				assert(fineMesh->Elements.size() == nx * ny * nz);
				if (nx == ny && ny == nz)
					assert(fineMesh->Faces.size() == 3 * n * n * (n + 1));
			}
			else if (meshCode.compare("stetra") == 0)
			{
				if (refinementStgy == H_CoarsStgy::StandardCoarsening || Utils::IsAlgebraic(args.Solver.SolverCode))
					fineMesh = new Cube_CartesianTetrahedralMesh(n);
				else
				{
					Mesh<3>* coarseMesh = new Cube_CartesianTetrahedralMesh(1);
					fineMesh = coarseMesh->RefineUntilNElements(6 * n * n * n, refinementStgy);
				}

				assert(fineMesh->Elements.size() == 6 * nx * ny * nz);
				if (nx == ny && ny == nz)
					assert(fineMesh->Faces.size() == 12 * n * n * n + 6 * n * n);
			}
			else if (meshCode.compare("tetra") == 0)
				Utils::FatalError("The in-house mesher does not build unstructured meshes. Use '-mesh stetra' or '-mesher gmsh'.");
			else
				Utils::FatalError("The requested mesh is not managed with this geometry.");
		}
#ifdef GMSH_ENABLED
		else if (mesher.compare("gmsh") == 0)
		{
			if (meshCode.compare("cart") == 0)
			{
				if (nx != ny || nx != nz)
					Utils::FatalError("Options -ny, -nz are not managed with this mesh");

				if (args.Solver.MG.H_CS == H_CoarsStgy::GMSHSplittingRefinement)
				{
					Mesh<3>* coarseMesh = new Cube_GMSHCartesianMesh(args.Solver.MG.CoarseN);
					fineMesh = coarseMesh->RefineUntilNElements(n * n * n, refinementStgy);

					assert(fineMesh->Elements.size() == n * n * n);
					assert(fineMesh->Faces.size() == 3 * n * n * (n + 1));
				}
				else
					fineMesh = new Cube_GMSHCartesianMesh(n);
			}
			else if (meshCode.compare("tetra") == 0)
			{
				if (nx != ny || nx != nz)
					Utils::FatalError("Options -ny, -nz are not managed with this mesh");

				if (Utils::IsRefinementStrategy(args.Solver.MG.H_CS))
				{
					Mesh<3>* coarseMesh = new Cube_GMSHTetrahedralMesh(args.Solver.MG.CoarseN);
					fineMesh = coarseMesh->RefineUntilNElements(6 * n * n * n, refinementStgy);

					assert(fineMesh->Elements.size() == 6 * n * n * n);
					assert(fineMesh->Faces.size() == 12 * n * n * n + 6 * n * n);
				}
				else
					fineMesh = new Cube_GMSHTetrahedralMesh(n);
			}
			else
				Utils::FatalError("The requested mesh is not managed with this geometry.");
		}
#endif // GMSH_ENABLED
	}
	//----------------------//
	//       GMSH file      //
	//----------------------//
	else
	{
		if (mesher.compare("inhouse") == 0)
			Utils::FatalError("The geometry is imported from a GMSH file, the mesher should be 'gmsh'. Use '-mesher gmsh' instead.");

#ifdef GMSH_ENABLED
		if (meshCode.compare("tetra") == 0)
		{
			string filePath = geoCode;
			if (nx != ny || nx != nz)
				Utils::FatalError("-ny, -ny not managed with this mesh");

			Mesh<3>* coarseMesh;
			if (Utils::IsRefinementStrategy(args.Solver.MG.H_CS))
			{
				if (refinementStgy == H_CoarsStgy::BeyRefinement)
					coarseMesh = new GMSHTetrahedralMesh(testCase, filePath, args.Solver.MG.CoarseN);
				else if (refinementStgy == H_CoarsStgy::GMSHSplittingRefinement)
					coarseMesh = new GMSHMesh<3>(testCase, filePath, args.Solver.MG.CoarseN);
				fineMesh = coarseMesh->RefineUntilNElements(6 * n * n * n, refinementStgy);
			}
			else
				fineMesh = new GMSHMesh<3>(testCase, filePath, n);
		}
		else
			Utils::FatalError("When the geometry is imported from a GMSH file, only unstructured tetrahedral meshing is allowed. Use '-mesh tetra' instead.");
#endif // GMSH_ENABLED
	}

	if (!Utils::IsRefinementStrategy(args.Solver.MG.H_CS) && fineMesh->CoarseMesh)
		fineMesh->DeleteCoarseMeshes();

	return fineMesh;
}
#endif // ENABLE_3D