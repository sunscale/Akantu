// Gmsh project created on Fri Dec 28 23:39:32 2018
Point(1) = {-0.5, -0.5, 0, 100.0};
Point(2) = {0.5, -0.5, 0, 100.0};
Point(3) = {0.5, 0.5, 0, 100.0};
Point(4) = {-0.5, 0.5, 0, 100.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {3, 4, 1, 2};

Plane Surface(1) = {1};
Physical Surface("plate") = {1};

Recombine Surface{1};