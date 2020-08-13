#pragma once
#include "1D/UniformMesh1D.h"

#include "2D/Square_CartesianMesh.h"
#include "2D/Square_CartesianPolygonalMesh.h"
#include "2D/Square_TriangularMesh.h"
#include "2D/Square_QuadrilateralMesh.h"
#include "2D/Square_QuadrilateralAsPolygonalMesh.h"

#include "3D/Cube_CartesianMesh.h"
#include "3D/Cube_CartesianTetrahedralMesh.h"

#ifdef GMSH_ENABLED
#include "2D/Square_GMSHCartesianMesh.h"
#include "2D/Square_GMSHTriangularMesh.h"
#include "2D/Square_GMSHUnstructTriangularMesh.h"
#include "2D/Square_GMSHQuadrilateralMesh.h"

#include "2D/Square4quadrants_GMSHCartesianMesh.h"
#include "2D/Square4quadrants_GMSHTriangularMesh.h"
#include "2D/Square4quadrants_GMSHUnstructTriangularMesh.h"
#include "2D/Square4quadrants_GMSHQuadrilateralMesh.h"

#include "3D/Cube_GMSHTetrahedralMesh.h"
#include "3D/Cube_GMSHCartesianMesh.h"
#endif