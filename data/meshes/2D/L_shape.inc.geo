L = 1; // size of the square
// h and h_center are defined in the files that include this one

// Corners
Point(1) = {0, 0, 0, h};
Point(2) = {0, L, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {L/2, L/2, 0, h_center}; // Center
Point(5) = {L/2, 0, 0, h};
Point(6) = {0, L/2, 0, h};
Point(7) = {L, L/2, 0, h};
Point(8) = {L/2, L, 0, h};

// Lines
Line(1) = {1, 5};
Line(2) = {7, 3};
Line(3) = {3, 8};
Line(4) = {8, 2};
Line(5) = {2, 6};
Line(6) = {6, 1};
Line(7) = {5, 4};
Line(8) = {4, 7};
Line(9) = {4, 8};
Line(10) = {6, 4};

// Define quadrants from boundaries
Curve Loop(1) = {1, 7, -10, 6};
Curve Loop(2) = {8, 2, 3, -9};
Curve Loop(3) = {10, 9, 4, 5};

// Define surfaces from quadrants
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// Define physical surface
Physical Surface("domain")  = {1, 2, 3};
