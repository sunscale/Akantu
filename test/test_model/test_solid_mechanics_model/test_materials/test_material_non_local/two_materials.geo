// Mesh size
h  = 0.25;    // Top cube

// Dimensions of top cube
Lx = 1;
Ly = 1;


// ------------------------------------------
// Geometry
// ------------------------------------------

// Base Cube
Point(101) = { 0.0, 0.0, 0.0, h}; // Bottom Face
Point(102) = { Lx/4, 0.0, 0.0, h};
Point(103) = { Lx,  0.0, 0.0, h}; // Bottom Face
Point(104) = { Lx,  Ly, 0.0,  h}; // Bottom Face
Point(105) = { Lx/4, Ly, 0.0, h};
Point(106) = { 0.0, Ly, 0.0,  h}; // Bottom Face

// Base Cube
Line(101) = {101,102}; // Bottom Face
Line(102) = {102,105}; // Bottom Face
Line(103) = {105,106}; // Bottom Face
Line(104) = {106,101}; // Bottom Face

Line(105) = {102,103}; // Bottom Face
Line(106) = {103,104}; // Bottom Face
Line(107) = {104,105}; // Bottom Face
// Line(108) = {105,102}; // Bottom Face
// 
// // Base Cube
// Line Loop(101) = {101:104};
// Line Loop(102) = {105:108};
// 
// Plane Surface(101) = {101};
// Plane Surface(102) = {102};
// Physical Surface(1) = {101};
// Physical Surface(2) = {102};

Line Loop(108) = {105, 106, 107, -102};
Plane Surface(209) = {108};
Line Loop(110) = {101, 102, 103, 104};
Plane Surface(210) = {110};

Physical Surface(1) = {209};
Physical Surface(2) = {210};

Transfinite Surface "*";
