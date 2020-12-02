L = 1; // size of the square
N = 10;
h = 1/(N+1);

// Big square's corners
Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};
// Small square's corners
Point(5) = {  L/4,   L/4, 0, h};
Point(6) = {3*L/4,   L/4, 0, h};
Point(7) = {3*L/4, 3*L/4, 0, h};
Point(8) = {  L/4, 3*L/4, 0, h};

// Big square's line
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
// Small square's line
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// Big square with the small square as a hole
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
// Small square
Curve Loop(2) = {5, 6, 7, 8};

Plane Surface(1) = {1};
Plane Surface(2) = {2};

// Define physical surfaces
Physical Surface("big")  = {1};
Physical Surface("small") = {2};
