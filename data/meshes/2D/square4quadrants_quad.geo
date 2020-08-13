N = 10;
h = 1/(N+1);

Include "square4quadrants.inc.geo";

// Combine triangles to make quadrilaterals
Recombine Surface "*";