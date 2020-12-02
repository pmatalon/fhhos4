L = 1; // size of the square
N = 10;
h = 1/(N+1);

// Corners
Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {0, L, 0, h};
Point(4) = {L, L, 0, h};
// Center
Point(5) = {L/2, L/2, 0, h};
// Middle edges
Point(6) = {L/2, 0, 0, h};
Point(7) = {L, L/2, 0, h};

// External lines
Line(1) = {1, 6};
Line(2) = {6, 2};
Line(3) = {2, 7};
Line(4) = {7, 4};
Line(5) = {4, 3};
Line(7) = {3, 1};
// Internal lines
Line(8) = {6, 5};
Line(9) = {5, 7};

// Squares' boundaries
Curve Loop(1) = {1, 8, 9, 4, 5, 7};
Curve Loop(2) = {2, 3, -9, -8};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

// Define physical surfaces
Physical Surface("big")  = {1};
Physical Surface("small") = {2};
