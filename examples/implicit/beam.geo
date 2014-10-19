h = 0.1;

Point(1) = {0,  0, 0, h};
Point(2) = {10, 0, 0, h};
Point(3) = {10, 1, 0, h};
Point(4) = {0,  1, 0, h};
Point(6) = {0, 0.5, 0, h};
Point(7) = {10, 0.5, 0, h};
Point(5) = {5, 0.5, 0, h};
Line(1) = {1, 2};
Line(21) = {2, 7};
Line(22) = {7, 3};
Line(3) = {3, 4};
Line(41) = {4, 6};
Line(42) = {6, 1};
Line Loop(5) = {3, 41, 42, 1, 21, 22};
Plane Surface(6) = {5};

Point{5} In Surface{6};

Physical Point ("traction") = {5};
Physical Point ("blocked")  = {6};
Physical Point ("roller")   = {7};
Physical Surface("bulk")    = {6};
Physical Line ("left")      = {41, 42};
Physical Line ("right")     = {21, 22};
