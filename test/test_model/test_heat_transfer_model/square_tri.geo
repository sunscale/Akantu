h = 0.01;

Point(1) = {0.0, 0.0, 0.0, h};
Point(2) = {1.0, 0.0, 0.0, h};
Point(3) = {0.0, 1.0, 0.0, h};
Point(4) = {1.0, 1.0, 0.0, h};

Line(3) = {4, 3};
Line(8) = {2, 1};
Line(9) = {3, 1};
Line(12) = {4, 2};

Line Loop(13) = {3, 9, -8, -12};
Plane Surface(14) = {13};

Transfinite Surface "*";
//Recombine Surface "*";
