h = 0.1;

Point(1) = {0.0, 0.0, 0.0, h};
Point(2) = {1.0, 0.0, 0.0, h};
Point(3) = {0.0, 1.0, 0.0, h};
Point(4) = {1.0, 1.0, 0.0, h};
Point(5) = {0.0, 0.0, 1.0, h};
Point(6) = {1.0, 0.0, 1.0, h};
Point(7) = {0.0, 1.0, 1.0, h};
Point(8) = {1.0, 1.0, 1.0, h};

Line(1) = {7, 8};
Line(2) = {8, 6};
Line(3) = {6, 5};
Line(4) = {5, 7};
Line(5) = {3, 7};
Line(6) = {5, 1};
Line(7) = {1, 3};
Line(8) = {2, 1};
Line(9) = {3, 4};
Line(10) = {2, 4};
Line(11) = {2, 6};
Line(12) = {8, 4};
Line Loop(13) = {3, 6, -8, 11};
Plane Surface(14) = {13};
Line Loop(15) = {4, -5, -7, -6};
Plane Surface(16) = {15};
Line Loop(17) = {1, 12, -9, 5};
Plane Surface(18) = {17};
Line Loop(19) = {1, 2, 3, 4};
Plane Surface(20) = {19};
Line Loop(21) = {12, -10, 11, -2};
Plane Surface(22) = {21};
Line Loop(23) = {9, -10, 8, 7};
Plane Surface(24) = {23};
Surface Loop(25) = {24, 18, 20, 22, 14, 16};
Volume(26) = {25};

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

Physical Volume (100) = {26};