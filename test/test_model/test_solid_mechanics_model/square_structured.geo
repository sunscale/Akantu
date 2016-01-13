// Mesh size
nb_nodes = 11.;

Lx = 1;
Ly = 1;

h = Lx / nb_nodes;

// ------------------------------------------
// Geometry
// ------------------------------------------

Point(1) = { 0.0, 0.0, 0.0, h};
Point(2) = { Lx,  0.0, 0.0, h};
Point(3) = { Lx,  Ly, 0.0,  h};
Point(4) = { 0.0, Ly, 0.0,  h};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Line Loop(1) = {1:4};
Plane Surface(2) = {1};

Transfinite Surface "*";
Recombine Surface "*";
