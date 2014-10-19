h=0.5;

Point(1) = {0,  0, 0, h};
Point(2) = {0, -1, 0, h};
Point(3) = {1, -1, 0, h};
Point(4) = {1,  0, 0, h};

Point(5) = {0,  0, 1, h};
Point(6) = {0, -1, 1, h};
Point(7) = {1, -1, 1, h};
Point(8) = {1,  0, 1, h};


Line(1) = {1, 5};
Line(2) = {5, 8};
Line(3) = {8, 4};
Line(4) = {4, 1};
Line(5) = {6, 2};
Line(6) = {2, 3};
Line(7) = {3, 7};
Line(8) = {7, 6};
Line(9) = {5, 6};
Line(10) = {1, 2};
Line(11) = {4, 3};
Line(12) = {8, 7};

Line Loop(26) = {4, 1, 2, 3};
Plane Surface(26) = {26};
Line Loop(28) = {9, -8, -12, -2};
Plane Surface(28) = {28};
Line Loop(30) = {-7, -11, -3, 12};
Plane Surface(30) = {30};
Line Loop(32) = {-10, -4, 11, -6};
Plane Surface(32) = {32};
Line Loop(34) = {10, -5, -9, -1};
Plane Surface(34) = {34};
Line Loop(36) = {8, 5, 6, 7};
Plane Surface(36) = {36};

Surface Loop(50) = {26, 28, 30, 32, 34, 36};
Volume(50) = {50};

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";

Physical Volume ("PhysVol") = {50};
Physical Surface ("Bottom") = {32}; // Bottom face
Physical Surface ("Top") = {28}; // Top face
Physical Surface ("Sides") = {26, 30, 34, 36}; // Side faces


