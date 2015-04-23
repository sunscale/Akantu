// Mesh size
h  = 0.05;    // Top cube

// Dimensions of top cube
Lx = 1;
Ly = 1;


// ------------------------------------------
// Geometry
// ------------------------------------------

// Base Cube
Point(101) = { 0.0, 0.0, 0.0, h}; // Bottom Face
Point(102) = { Lx,  0.0, 0.0, h}; // Bottom Face
Point(103) = { Lx,  Ly, 0.0,  h}; // Bottom Face
Point(104) = { 0.0, Ly, 0.0,  h}; // Bottom Face

// Base Cube
Line(101) = {101,102}; // Bottom Face
Line(102) = {102,103}; // Bottom Face
Line(103) = {103,104}; // Bottom Face
Line(104) = {104,101}; // Bottom Face

// Base Cube
Line Loop(101) = {101:104};

Plane Surface(101) = {101};

Physical Line("Fixed") = {101, 103, 104};
Physical Line("Traction") = {102};
Physical Surface("Interior") = {101};

//Transfinite Surface "*";
//Recombine Surface "*";
