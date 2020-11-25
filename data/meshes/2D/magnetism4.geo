L = 1;
N = 8;
h = 1/(N+1);

Point(1)   = {0     , 0, 0, h}; // center
Point(2)   = {     L, 0, 0, h};
Point(6)   = { 0.4*L, 0, 0, h};
Point(39)  = {-0.4*L, 0, 0, h};
Point(100) = {    -L, 0, 0, h};

// Exterior half-circles
Circle(145) = {2, 1, 100};
Circle(146) = {100, 1, 2};

// Interior half-circles
Circle(151) = {6, 1, 39};
Circle(152) = {39, 1, 6};

Line Loop(149) = {145, 146, -151, -152};
Plane Surface(149) = {149};

// Interior circle
Line Loop(156) = {151, 152};
Physical Line("innerCircle") = {156};
Plane Surface(156) = {156};

Physical Line("externalBoundary") = {145, 146};

Physical Surface("Exterior") = {149};
Physical Surface("Interior") = {156};
