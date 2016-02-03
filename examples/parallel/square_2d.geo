// Element size
h = 0.05;

// Dimension of square
L = 1;

// ------------------------------------------
// Geometry
// ------------------------------------------

// Points
Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, L, 0, h};
Point(4) = {0, L, 0, h};

// Lines
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};

// Geometric and Physical Surface
Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Physical Surface(7) = {6};
