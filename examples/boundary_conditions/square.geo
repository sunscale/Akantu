// Mesh size
h  = 0.005;

// Dimensions of the square
Lx = 0.01;
Ly = 0.01;

// ------------------------------------------
// Geometry
// ------------------------------------------
Point(1) = { 0.0, 0.0, 0.0, h};
Point(2) = { Lx,  0.0, 0.0, h};
Point(3) = { Lx,  Ly, 0.0,  h};
Point(4) = { 0.0, Ly, 0.0,  h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1:4};

Plane Surface(1) = {1};

Physical Surface(1) = {1};
Physical Line("Fixed_y")  = {1};
Physical Line("Fixed_x")  = {4};
Physical Line("Traction") = {2};
Physical Line("Free")     = {3};
