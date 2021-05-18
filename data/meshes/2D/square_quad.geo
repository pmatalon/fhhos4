N = 10;
h = 1/(N+1);
h_ref = h;

Include "square.inc.geo";

// Combine triangles to make quadrilaterals
Recombine Surface "*";