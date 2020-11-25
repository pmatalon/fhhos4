N=4;
h=1/N;
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {0.8, 0, 0, h};
Point(6) = {0.4, 0, 0, h};

Point(37) = {-0.8, 0, 0, h};
Point(39) = {-0.4, 0, 0, h};
Point(100) = {-1, 0, 0, h};
Circle(10) = {3, 1, 37};
Circle(15) = {37, 1, 3};
Circle(145) = {2, 1, 100};
Circle(146) = {100, 1, 2};
Circle(151) = {6, 1, 39};
Circle(152) = {39, 1, 6};
Line Loop(149) = {145, 146, -10, -15};
Plane Surface(149) = {149};
Line Loop(155) = {15, 10, -152, -151};
Plane Surface(155) = {155};
Line Loop(156) = {151, 152};
Plane Surface(156) = {156};
Physical Line(1) = {145, 146};
Physical Surface("Exterior") = {149};
Physical Surface("Middle") = {155};
Physical Surface("Interior") = {156};
