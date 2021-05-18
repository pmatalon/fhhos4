L = 1;
// h is defined in the files that include this one

// Corners
Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h_ref};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};

// External lines
Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 1}; // left

// Boundaries
Curve Loop(1) = {1, 2, 3, 4};

// Surface
Plane Surface(1) = {1};

Physical Line("bottomBoundary") = {1};
Physical Line("rightBoundary")  = {2};
Physical Line("topBoundary")    = {3};
Physical Line("leftBoundary")   = {4};

Physical Surface("domain") = {1};