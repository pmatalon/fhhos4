// Gmsh project created on Fri Dec 13 13:50:13 2019
//SetFactory("OpenCASCADE");
L = 1;
N = 2;
h = 1/N;

// Corners
Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {0, L, 0, h};
Point(4) = {L, L, 0, h};
// Center
Point(5) = {L/2, L/2, 0, h};
// Middle edges
Point(6) = {L/2, 0, 0, h};
Point(7) = {0, L/2, 0, h};
Point(8) = {L, L/2, 0, h};
Point(9) = {L/2, L, 0, h};

// External lines
Line(1) = {1, 6};
Line(2) = {6, 2};
Line(3) = {2, 8};
Line(4) = {8, 4};
Line(5) = {4, 9};
Line(6) = {9, 3};
Line(7) = {3, 7};
Line(8) = {7, 1};
// Internal lines
Line(9) = {6, 5};
Line(10) = {5, 8};
Line(11) = {5, 9};
Line(12) = {7, 5};

// Quadrants boundaries
Curve Loop(1) = {1, 9, -12, 8};
Curve Loop(2) = {2, 3, -10, -9};
Curve Loop(3) = {10, 4, 5, -11};
Curve Loop(4) = {12, 11, 6, 7};

// Quadrant surfaces from their boundaries
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
