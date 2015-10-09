
// Mesh size
h  = .025;
// Dimensions of the square
Lx = 2;
Ly = .2;
// corners
Point(1) = { 0.0, 0.0, 0.0, h}; // Bottom Face
Point(2) = { Lx,  0.0, 0.0, h}; // Bottom Face
Point(3) = { Lx,  Ly, 0.0,  h}; // Bottom Face
Point(4) = { 0.0, Ly, 0.0,  h}; // Bottom Face
// edges
Line(101) = {1,2}; // Bottom Face
Line(102) = {2,3}; // Bottom Face
Line(103) = {3,4}; // Bottom Face
Line(104) = {4,1}; // Bottom Face
// square perimeter
Line Loop(101) = {101:104};
// square surface
Plane Surface(101) = {101};
// physical surface
Physical Surface("Bulk") = {101};
