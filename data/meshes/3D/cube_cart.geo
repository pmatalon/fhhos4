N = 2;
h = 1/(N+1);

Include "cube.inc.geo";

// Structured tetrahedral
Transfinite Line {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12} = N+1;
Transfinite Surface "*";
Transfinite Volume "*";

// Combine tetra to create cubes
Recombine Surface "*";
Recombine Volume "*";
