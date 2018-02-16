// Mesh size
h  = 0.05;

h1 = h;
h2 = h;

// Dimensions of the bar
Lx = 10;
Ly = 1;

// ------------------------------------------
// Geometry
// ------------------------------------------

Point(101) = { 0.0, -Ly/2, 0.0, h1};
Point(102) = { Lx,  -Ly/2, 0.0, h2};

Point(103) = { Lx,  0., 0.0,  h2};
Point(104) = { Lx,  Ly/2., 0.0,  h2};

Point(105) = { 0.0, Ly/2., 0.0,  h1};
Point(106) = { 0.0, 0., 0.0,  h1};

Line(101) = {101,102};
Line(102) = {102,103};
Line(103) = {103,104};
Line(104) = {104,105};
Line(105) = {105,106};
Line(106) = {106,101};
Line(107) = {106,103};


Line Loop(108) = {101, 102, -107, 106};
Plane Surface(109) = {108};
Line Loop(110) = {103, 104, 105, 107};
Plane Surface(111) = {110};
Physical Surface(112) = {109, 111};

Transfinite Surface "*";
Recombine Surface "*";
Physical Surface(113) = {111, 109};

Physical Line("XBlocked") = {103, 102};
Physical Line("YBlocked") = {104, 101};
