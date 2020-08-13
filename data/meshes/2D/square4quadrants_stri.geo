N = 10;
h = 1/(N+1);

Include "square4quadrants.inc.geo";

// Structured triangular
Transfinite Line {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12} = N+1;
Transfinite Surface "*";