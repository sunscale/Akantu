// Mesh size
h  = .05;

h1 = h;
h2 = h;

// Dimensions of the bar
Lx = 1.;
Ly = 1.;

// ------------------------------------------
// Geometry
// ------------------------------------------

Point(101) = { 0.0, 0.0, 0.0, h1};
Point(102) = { Lx,  0.0, 0.0, h2};

Point(103) = { Lx,  Ly, 0.0,  h2};
Point(104) = { 0.,  Ly, 0.0,  h2};
//+
Line(1) = {101, 102};
//+
Line(2) = {102, 103};
//+
Line(3) = {103, 104};
//+
Line(4) = {104, 101};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Surface("bulk") = {1};
