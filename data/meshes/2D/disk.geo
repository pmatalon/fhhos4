L = 1;
N = 8;
h = 1/(N+1);

Point(1)   = {0     , 0, 0, h}; // center
Point(2)   = {     L, 0, 0, h};
Point(6)   = { 0.4*L, 0, 0, h};
Point(39)  = {-0.4*L, 0, 0, h};
Point(100) = {    -L, 0, 0, h};

// Half-circles
Circle(145) = {2, 1, 100};
Circle(146) = {100, 1, 2};

Line Loop(149) = {145, 146};
Plane Surface(149) = {149};

Physical Line("externalBoundary") = {145, 146};

Physical Surface("Domain") = {149};
