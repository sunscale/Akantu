lc = 1;
h = 1;
Point(1) = {lc, lc, lc, h};
Point(2) = {-lc, lc, lc, h};
Point(3) = {-lc, -lc, lc, h};
Point(4) = {lc, -lc, lc, h};
Point(5) = {lc, lc, -lc, h};
Point(6) = {-lc, lc, -lc, h};
Point(7) = {-lc, -lc, -lc, h};
Point(8) = {lc, -lc, -lc, h};
Point(9) = {0, lc, lc, h};
Point(10) = {0, -lc, lc, h};
Point(11) = {0, -lc, -lc, h};
Point(12) = {0, lc, -lc, h};
Line(1) = {2, 3};
Line(2) = {3, 10};
Line(3) = {10, 9};
Line(4) = {9, 2};
Line(5) = {10, 4};
Line(6) = {4, 1};
Line(7) = {1, 9};
Line(8) = {10, 11};
Line(9) = {4, 8};
Line(10) = {3, 7};
Line(11) = {7, 11};
Line(12) = {11, 8};
Line(13) = {8, 5};
Line(14) = {5, 1};
Line(15) = {12, 9};
Line(16) = {6, 2};
Line(17) = {7, 6};
Line(18) = {6, 12};
Line(19) = {12, 5};
Line(20) = {11, 12};
Line Loop(21) = {8, 20, 15, -3};
Plane Surface(22) = {21};
Line Loop(23) = {9, 13, 14, -6};
Plane Surface(24) = {23};
Line Loop(25) = {10, 17, 16, 1};
Plane Surface(26) = {25};
Line Loop(27) = {8, -11, -10, 2};
Plane Surface(28) = {27};
Line Loop(29) = {9, -12, -8, 5};
Plane Surface(30) = {29};
Line Loop(31) = {7, -15, 19, 14};
Plane Surface(32) = {31};
Line Loop(33) = {4, -16, 18, 15};
Plane Surface(34) = {33};
Line Loop(35) = {5, 6, 7, -3};
Plane Surface(36) = {35};
Line Loop(37) = {2, 3, 4, 1};
Plane Surface(38) = {37};
Line Loop(39) = {18, -20, -11, 17};
Plane Surface(40) = {39};
Line Loop(41) = {19, -13, -12, 20};
Plane Surface(42) = {41};
Surface Loop(43) = {34, 38, 28, 40, 26, 22};
Volume(44) = {43};
Surface Loop(45) = {42, 32, 36, 30, 24, 22};
Volume(46) = {45};

Physical Volume(47) = {44};
Physical Volume(48) = {46};