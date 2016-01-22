cl1 = 4;
cl2 = 1;
offx = 0.;
offz = 0.25;
h = -1.3;
Point(1) = {-0.5, -0.5, -0.5, cl1};
Point(2) = {-0.5, 0.5, -0.5, cl1};
Point(3) = {1, 0.5, -0.5, cl1};
Point(4) = {1, -0.5, -0.5, cl1};
Point(5) = {-0.15+offx, 1.9+h, 0.4 + offz, cl2};
Point(6) = {-0.15+offx, 1.9+h, -0.4 + offz, cl2};
Point(7) = {0.65+offx, 1.9+h, -0.4 + offz, cl2};
Point(8) = {0.65+offx, 1.9+h, 0.4 + offz, cl2};
Point(11) = {-0.5, -0.5, 1, cl1};
Point(12) = {-0.5, 0.5, 1, cl1};
Point(13) = {1, 0.5, 1, cl1};
Point(14) = {1, -0.5, 1, cl1};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line(11) = {11, 14};
Line(12) = {14, 13};
Line(13) = {13, 12};
Line(14) = {12, 11};
Line(15) = {14, 4};
Line(16) = {13, 3};
Line(17) = {11, 1};
Line(18) = {12, 2};
Line(19) = {8, 7};
Line(20) = {7, 6};
Line(21) = {6, 5};
Line(22) = {5, 8};
Line Loop(23) = {22, 19, 20, 21};
Plane Surface(24) = {23};
Extrude {0, 1.5, 0} {
  Surface{24};
}
Coherence;

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

Surface Loop(50) = {36, 26, 32, 34, 28, 30};
Volume(50) = {50};


Physical Surface("Top") = {46};
Physical Surface("Bottom") = {34};
Physical Surface("Contact") = {36};
Physical Volume("Bottom_body") = {50};
Physical Volume("Top_body") = {1};
