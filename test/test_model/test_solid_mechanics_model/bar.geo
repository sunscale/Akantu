// Mesh size
h  = 0.2;    // Top cube

// Dimensions of top cube
Lx = 10;
Ly = 1;


// ------------------------------------------
// Geometry
// ------------------------------------------

// Base Cube
Point(101) = { 0.0, 0.0, 0.0, h}; //Left  bottom corner
Point(102) = { Lx,  0.0, 0.0, h}; //Right bottom corner
Point(103) = { Lx,  Ly, 0.0,  h}; //Right top    corner
Point(104) = { 0.0, Ly, 0.0,  h}; //Left  top    corner

// Base Cube
Line(101) = {101,102}; // Bottom Face
Line(102) = {102,103}; // Right  Face
Line(103) = {103,104}; // Top    Face
Line(104) = {104,101}; // Left   Face

// Base Cube
Line Loop(101) = {101:104};
Plane Surface(101) = {101};
Physical Surface(1) = {101};

Physical Line("Right")  = {102};
Physical Line("Left")   = {104};
Physical Line("Top")    = {103};
Physical Line("Bottom") = {101};

Extrude {0, 0, 1} {
  Surface{101};
}
