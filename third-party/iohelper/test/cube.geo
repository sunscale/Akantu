h = 1.;


Point(1) = {0.0, 0.0, 0.0, h};
Point(2) = {1.0, 0.0, 0.0, h};
Point(3) = {0.0, 1.0, 0.0, h};
Point(4) = {1.0, 1.0, 0.0, h};
Point(5) = {0.0, 0.0, 1.0, h};
Point(6) = {1.0, 0.0, 1.0, h};
Point(7) = {0.0, 1.0, 1.0, h};
Point(8) = {1.0, 1.0, 1.0, h};

Line(1) = {7, 5};
Line(2) = {5, 6};
Line(3) = {6, 8};
Line(4) = {6, 2};
Line(5) = {5, 1};
Line(6) = {7, 3};
Line(7) = {3, 1};
Line(8) = {7, 8};
Line(9) = {8, 4};
Line(10) = {4, 2};
Line(11) = {2, 1};
Line(12) = {3, 4};
Line Loop(13) = {6, 7, -5, -1};
Plane Surface(14) = {13};
Line Loop(15) = {8, -3, -2, -1};
Plane Surface(16) = {15};
Line Loop(17) = {4, -10, -9, -3};
Plane Surface(18) = {17};
Line Loop(19) = {12, 10, 11, -7};
Plane Surface(20) = {19};
Line Loop(21) = {11, -5, 2, 4};
Plane Surface(22) = {21};
Line Loop(23) = {6, 12, -9, -8};
Plane Surface(24) = {23};
Surface Loop(25) = {14, 24, 20, 18, 22, 16};

Volume(26) = {25};

Physical Volume(27) = {26};
