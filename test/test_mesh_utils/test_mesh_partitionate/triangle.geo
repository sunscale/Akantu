// Mesh size
h  = .1;    // Top cube

// Dimensions of top cube
Lx = 1;
Ly = 1;

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

// Base Cube
Plane Surface(101) = {101};

Physical Surface(6) = {101};

