N = 10;
h = 1/(N+1);

// Room corners
L = 1;
Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};

// Wall corners
T = 0.02; // thickness
G = 0.75; // length
Y = 0.3; // Y location
Point(5) = {0, Y,   0, h};
Point(6) = {G, Y,   0, h};
Point(7) = {G, Y+T, 0, h};
Point(8) = {0, Y+T, 0, h};

// External lines
Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 8}; // left-top
Line(5) = {8, 5}; // wall-left
Line(6) = {5, 1}; // left-bottom
// Internal lines
Line(7) = {5, 6}; // wall-bottom
Line(8) = {6, 7}; // wall-right
Line(9) = {7, 8}; // wall-top

// Air
Curve Loop(1) = {1, 2, 3, 4, -9, -8, -7, 6};
// Wall
Curve Loop(2) = {7, 8, 9, 5};

// Surfaces
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Physical Line("bottomBoundary") = {1};
Physical Line("rightBoundary")  = {2};
Physical Line("topBoundary")    = {3};
Physical Line("leftBoundary")   = {4, 5, 6};

Physical Surface("air") = {1};
Physical Surface("wall") = {2};