// CUBE ON CUBE (2 BODIES) 3D:

// REMARKS:
//   Physical Surfaces are defined so that load_mesh_msh can load
// directly the surface element structure. When creating Plane Surface,
// the surface normal has to point to the inside of the body. Otherwise the
// solvecontact3d algorithm won't work. The Physical Surface number 
// corresponds to the face number. The face numbering has to start with 1!
//   Physical Volume defines the material of each body. This can be read
// by the load_mesh_msh function in adlib.

// Mesh size
h  = 0.02;    // Top cube
hy = 0.5;    // Top cube

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

// Base Cube
Plane Surface(101) = {101};
Physical Line("Edges") = {101, 102, 103, 104};
Physical Surface(6) = {101};
