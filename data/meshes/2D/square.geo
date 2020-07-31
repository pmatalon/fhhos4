L = 1;
N = 10;
h = 1/N;

// Corners
Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};

// External lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Boundaries
Curve Loop(1) = {1, 2, 3, 4};

// Surface from boundaries
Plane Surface(1) = {1};
