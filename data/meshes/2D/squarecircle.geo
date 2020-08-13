L = 1;
N = 10;
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
Physical Line(1) = {1, 2, 3, 4};

// Circle points
Point(5) = {L/2, L/2, 0, h}; // Center
Point(6) = {L/4, L/2, 0, h}; // Left
Point(7) = {3*L/4, L/2, 0, h}; // Right
//Point(8) = {L/2, 3*L/4, 0, h}; // Top
//Point(9) = {L/2, L/4, 0, h}; // Bottom

// 4 cicle arcs
//Circle(10) = {7, 5, 8};
//Circle(11) = {8, 5, 6};
//Circle(12) = {6, 5, 9};
//Circle(13) = {9, 5, 7};
Circle(10) = {6, 5, 7};
Circle(11) = {7, 5, 6};

// Circle
//Line Loop(6) = {10, 11, 12, 13};
Line Loop(6) = {10, 11};
Physical Line(2) = {6};

Transfinite Line {1, 2, 3, 4} = 4;

Plane Surface(14) = {5, 6};
Physical Surface(1) = {14};

Plane Surface(15) = {6};
Physical Surface(2) = {15};