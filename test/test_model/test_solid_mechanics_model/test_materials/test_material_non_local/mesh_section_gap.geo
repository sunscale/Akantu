// Mesh size
h  = 0.05;

// Dimensions of thin part
Lx = 1;
Ly = 0.2;


// ------------------------------------------
// Geometry
// ------------------------------------------

// Base Cube
Point(101) = { 0.0, 0.0, 0.0, h}; // Bottom Face
Point(102) = { Lx,  0.0, 0.0, h}; // Bottom Face
Point(103) = { Lx,  2*Ly, 0.0,  h}; 
Point(104) = { Lx/2, 2*Ly, 0.0,  h};
Point(105) = { Lx/2, Ly, 0.0,  h};
Point(106) = { 0.0, Ly, 0.0,  h};

// Base Cube
Line(101) = {101,102}; // Bottom Face
Line(102) = {102,103};
Line(103) = {103,104};
Line(104) = {104,105};
Line(105) = {105,106};
Line(106) = {106,101};

// Base Cube
Line Loop(101) = {101:106};

Plane Surface(101) = {101};

Physical Line("Fixed") = {101, 106};
Physical Line("Traction") = {102};
Physical Surface("Interior") = {101};

//Transfinite Surface "*";
//Recombine Surface "*";
