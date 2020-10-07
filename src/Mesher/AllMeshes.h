#pragma once
#include "InHouse/UniformMesh1D.h"

#include "InHouse/Square_CartesianMesh.h"
#include "InHouse/Square_CartesianPolygonalMesh.h"
#include "InHouse/Square_TriangularMesh.h"
#include "InHouse/Square_QuadrilateralMesh.h"
#include "InHouse/Square_QuadrilateralAsPolygonalMesh.h"

#include "InHouse/Cube_CartesianMesh.h"
#include "InHouse/Cube_CartesianTetrahedralMesh.h"

#ifdef GMSH_ENABLED
#include "GMSH/Square_GMSHCartesianMesh.h"
#include "GMSH/Square_GMSHTriangularMesh.h"
#include "GMSH/Square_GMSHUnstructTriangularMesh.h"
#include "GMSH/Square_GMSHQuadrilateralMesh.h"

#include "GMSH/Square4quadrants_GMSHCartesianMesh.h"
#include "GMSH/Square4quadrants_GMSHTriangularMesh.h"
#include "GMSH/Square4quadrants_GMSHUnstructTriangularMesh.h"
#include "GMSH/Square4quadrants_GMSHQuadrilateralMesh.h"

#include "GMSH/Cube_GMSHTetrahedralMesh.h"
#include "GMSH/Cube_GMSHCartesianMesh.h"
#endif