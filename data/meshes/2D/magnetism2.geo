N=4;
h=1/N;
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {0.8, 0, 0, h};
Point(6) = {0.4, 0, 0, h};
// Points B et B'
y=Sqrt(0.8^2-0.2^2);
Point(21) = {0.2, y, 0, h};
Point(31) = {-0.2, y, 0, h};
// Points D et D'
y=Sqrt(0.6^2-0.2^2);
Point(23) = {0.2, y, 0, h};
Point(33) = {-0.2, y, 0, h};
// Points E et E'
y=Sqrt(0.5^2-0.2^2);
Point(24) = {0.2, y, 0, h};
Point(34) = {-0.2, y, 0, h};

Point(37) = {-0.8, 0, 0, h};
Point(39) = {-0.4, 0, 0, h};
Point(100) = {-1, 0, 0, h};
Line(6) = {23, 21};
Line(9) = {31, 33};
Circle(10) = {34, 1, 24};
Line(11) = {34, 33};
Line(12) = {24, 23};
Circle(13) = {21, 1, 3};
Circle(14) = {31, 1, 37};
Circle(15) = {37, 1, 3};
Circle(145) = {2, 1, 100};
Circle(146) = {100, 1, 2};
Circle(151) = {6, 1, 39};
Circle(152) = {39, 1, 6};
Line Loop(149) = {145, 146, 9, -11, 10, 12, 6, 13, -15, -14};
Plane Surface(149) = {149};
Line Loop(155) = {14, 15, -13, -6, -12, -10, 11, -9, -152, -151};
Plane Surface(155) = {155};
Line Loop(156) = {151, 152};
Plane Surface(156) = {156};
Physical Line(1) = {145, 146};
Physical Surface("Exterior") = {149};
Physical Surface("Middle") = {155};
Physical Surface("Interior") = {156};
