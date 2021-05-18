N = 10;
h = 1/(N+1);
h_ref = h;

Include "square.inc.geo";

// Structured triangular
Transfinite Line {1, 2, 3, 4} = N+1;
Transfinite Surface "*";
