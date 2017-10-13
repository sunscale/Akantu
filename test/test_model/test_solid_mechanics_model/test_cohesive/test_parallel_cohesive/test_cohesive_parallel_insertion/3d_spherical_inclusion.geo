dx = 1;
w = 2;
radius = 0.2;

Point(1) = {0,0,0,dx};
Point(2) = {0,.3,0,dx};
Point(3) = {0,-.3,0,dx};
Point(4) = {w,0,0,dx};
Point(5) = {w,.3,0,dx};
Point(6) = {w,-.3,0,dx};
Point(7) = {0,0,w,dx};
Point(8) = {0,.3,w,dx};
Point(9) = {0,-.3,w,dx};
Point(10) = {w,0,w,dx};
Point(11) = {w,.3,w,dx};
Point(12) = {w,-.3,w,dx};
Point(13) = {0.5*w,0,0.5*w,dx};
Point(14) = {0.5*w+radius,0,0.5*w+radius,dx};
Point(15) = {0.5*w-radius,0,0.5*w+radius,dx};
Point(16) = {0.5*w+radius,0,0.5*w-radius,dx};
Point(17) = {0.5*w-radius,0,0.5*w-radius,dx};
Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 4};
Line(4) = {1, 4};
Line(5) = {1, 3};
Line(6) = {6, 4};
Line(7) = {3, 6};
Line(8) = {8, 11};
Line(9) = {7, 10};
Line(10) = {9, 12};
Line(11) = {8, 7};
Line(12) = {11, 10};
Line(13) = {10, 12};
Line(14) = {7, 9};
Line(15) = {2, 8};
Line(16) = {1, 7};
Line(17) = {3, 9};
Line(18) = {5, 11};
Line(19) = {4, 10};
Line(20) = {6, 12};
Line Loop(21) = {18, 12, -19, -3};
Plane Surface(22) = {21};
Line Loop(23) = {19, 13, -20, 6};
Plane Surface(24) = {23};
Line Loop(25) = {11, -16, 1, 15};
Plane Surface(26) = {25};
Line Loop(27) = {14, -17, -5, 16};
Plane Surface(28) = {27};
Line Loop(29) = {8, 12, -9, -11};
Plane Surface(30) = {29};
Line Loop(31) = {9, 13, -10, -14};
Plane Surface(32) = {31};
Line Loop(33) = {2, 3, -4, 1};
Plane Surface(34) = {33};
Line Loop(35) = {7, 6, -4, 5};
Plane Surface(36) = {35};
Line Loop(37) = {18, -8, -15, 2};
Plane Surface(38) = {37};
Line Loop(39) = {7, 20, -10, -17};
Plane Surface(40) = {39};
Circle(41) = {16, 13, 14};
Circle(42) = {14, 13, 15};
Circle(43) = {15, 13, 17};
Circle(44) = {17, 13, 16};
Translate {0.5, 0, 0.5} {
  Duplicata { Line{44, 43, 42, 41}; }
}
Translate {-0.5, 0, 0.5} {
  Duplicata { Line{44, 43, 42, 41}; }
}
Translate {0.5, 0, -0.5} {
  Duplicata { Line{44, 43, 42, 41}; }
}
Translate {-0.5, 0, -0.5} {
  Duplicata { Line{44, 43, 42, 41}; }
}
Line Loop(61) = {9, -19, -4, 16};
Line Loop(62) = {52, 51, 50, 49};
Line Loop(63) = {57, 60, 59, 58};
Line Loop(64) = {44, 41, 42, 43};
Line Loop(65) = {55, 54, 53, 56};
Line Loop(66) = {48, 47, 46, 45};
Plane Surface(67) = {61, 62, 63, 64, 65, 66};
Plane Surface(68) = {62};
Plane Surface(69) = {63};
Plane Surface(70) = {64};
Plane Surface(71) = {66};
Plane Surface(72) = {65};
Surface Loop(73) = {26, 30, 38, 22, 34, 67, 68, 69, 70, 72, 71};
Volume(74) = {73};
Surface Loop(75) = {32, 24, 40, 36, 28, 67, 68, 69, 70, 72, 71};
Volume(76) = {75};
Physical Surface("interface") = {67};
Physical Surface("coh1") = {68};
Physical Surface("coh2") = {69};
Physical Surface("coh3") = {70};
Physical Surface("coh4") = {71};
Physical Surface("coh5") = {72};
Physical Volume("bulk") = {74, 76};
