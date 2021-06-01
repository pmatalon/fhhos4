N = 5;
h = 1/(N+1);

Point(1) = {-10, 10, 0, 10*h};
Point(2) = {-10, -10, 0, 10*h};
Point(3) = {10, -10, 0, 10*h};
Point(4) = {10, 10, 0, 10*h};

Point(10) = {-8, 1, 0, h};
Point(11) = {-8, -1, 0, h};
Point(12) = {8, -1, 0, h};
Point(13) = {8, 1, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {10, 11};
Line(6) = {11, 12};
Line(7) = {12, 13};
Line(8) = {13, 10};

Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,8};

Plane Surface(1) = {1,2};
Plane Surface(2) = {2};

Transfinite Line {6,8} = 10*N;
Transfinite Line {5,7} = N+2;

Transfinite Surface{2} = {10,11,12,13};

Recombine Surface{2};


Physical Surface("exterior") = {1};
Physical Surface("interior") = {2};
