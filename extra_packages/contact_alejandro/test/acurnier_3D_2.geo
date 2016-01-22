cl1 = 4;
off = 0.25;
h = -1.3;
Point(1) = {-0.5, -0.5, -0.5, cl1};
Point(2) = {-0.5, 0.5, -0.5, cl1};
Point(3) = {1, 0.5, -0.5, cl1};
Point(4) = {1, -0.5, -0.5, cl1};
Point(5) = {0.75, 2.8+h, -0.5 + off, cl1};
Point(6) = {0.75, 1.9+h, -0.5 + off, cl1};
Point(7) = {-0.25, 1.9+h, -0.5 + off, cl1};
Point(11) = {-0.5, -0.5, 1, cl1};
Point(12) = {-0.5, 0.5, 1, cl1};
Point(13) = {1, 0.5, 1, cl1};
Point(14) = {1, -0.5, 1, cl1};
Point(26) = {0.75, 1.9+h, 0.5 + off, cl1};
Point(30) = {-0.25, 2.8+h, -0.5 + off, cl1};
Point(31) = {0.75, 2.8+h, 0.5 + off, cl1};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 5};
Line(11) = {11, 14};
Line(12) = {14, 13};
Line(13) = {13, 12};
Line(14) = {12, 11};
Line(15) = {14, 4};
Line(16) = {13, 3};
Line(17) = {11, 1};
Line(18) = {12, 2};
Line(19) = {7, 26};
Line(21) = {26, 6};
Line(23) = {26, 5};

Line Loop(26) = {-12, 15, 2, -16};
Plane Surface(26) = {26};
Line Loop(28) = {-1, -4, -3, -2};
Plane Surface(28) = {28};
Line Loop(30) = {-17, -14, 18, 4};
Plane Surface(30) = {30};
Line Loop(32) = {11, 12, 13, 14};
Plane Surface(32) = {32};
Line Loop(34) = {-15,1, 17, -11};
Plane Surface(34) = {34};
Line Loop(36) = {16, 3, -18, -13};
Plane Surface(36) = {36};

Line Loop(38) = {-5, -23, 21};
Plane Surface(38) = {38};
Line Loop(44) = {7, 6, 5};
Plane Surface(44) = {44};
Line Loop(45) = {23, 19, -7};
Plane Surface(45) = {45};
Line Loop(46) = {-21, -19, -6};
Plane Surface(46) = {46};

Line(51) = {31, 5};
Line(52) = {5, 30};
Line(53) = {30, 31};
Line(54) = {31, 26};
Line(55) = {30, 7};
Line(67) = {30, 26};

Line Loop(56) = {51, 52, 53};
Plane Surface(57) = {56};
Line Loop(58) = {54, 23, -51};
Plane Surface(59) = {58};


Line Loop(62) = {-7, -55, -52};
Plane Surface(63) = {62};

Line Loop(68) = {-54, -53, 67};
Plane Surface(69) = {68};
Line Loop(70) = {-67, 55, 19};
Plane Surface(71) = {70};
Line Loop(72) = {52, 67, 23};
Plane Surface(73) = {72};


Surface Loop(50) = {36, 26, 32, 34, 28, 30};
Volume(50) = {50};

Surface Loop(79) = {38, 44, 46, 63, 71, 69, 59, 57};
Volume(80) = {79};

Physical Surface("Top") = {57};
Physical Surface("Bottom") = {34};
Physical Surface("Contact") = {36};
Physical Volume("Top_body") = {80};
Physical Volume("Bottom_body") = {50};

