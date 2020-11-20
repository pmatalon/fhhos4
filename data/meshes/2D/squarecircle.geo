L = 1;
N = 8;
h = 1/(N+1);

// Square
Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};

// Square lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};

// Boundary
Physical Line("bottomBoundary") = {1};
Physical Line("rightBoundary")  = {2};
Physical Line("topBoundary")    = {3};
Physical Line("leftBoundary")   = {4};

// Circle points
Point(5) = {L/2, L/2, 0, h}; // Center
Point(6) = {L/4, L/2, 0, h}; // Left
Point(7) = {3*L/4, L/2, 0, h}; // Right

// 2 circle arcs
Circle(10) = {6, 5, 7};
Circle(11) = {7, 5, 6};

// Circle
Line Loop(6) = {10, 11};
Physical Line("circle") = {6};

Plane Surface(14) = {5, 6};
Physical Surface("square") = {14};

Plane Surface(15) = {6};
Physical Surface("disk") = {15};