// Element size
h = 0.5;

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

Point(5) = {1.5*L, 0, 0, h};
Point(6) = {2.5*L, 0, 0, h};
Point(7) = {2.5*L, L, 0, h};
Point(8) = {1.5*L, L, 0, h};

// Lines
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};

Line(5) = {8, 5};
Line(6) = {5, 6};
Line(7) = {6, 7};
Line(8) = {7, 8};

// Geometric and Physical Surface
Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
