cl1 = 0.02;
cl2 = 0.2;
cl3 = 0.4;
Dy = 0.099;
Dz = 1;
Point(1) = {0, 0.1-Dy, 0, cl1};
Point(2) = {0.5, 0.6-Dy, 0, cl2};
Point(3) = {-0.5, 0.6-Dy, 0, cl2};
Point(4) = {0, 0.6-Dy, 0, cl2};
Point(5) = {0, 0.6-Dy, -0.5, cl2};
Point(6) = {0, 0.6-Dy, 0.5, cl2};

Point(7) = {0, 0, 0, cl3};
Point(8) = {0.5, 0, 0, cl3};
Point(9) = {-0.5, 0, 0, cl3};
Point(10) = {0, 0, -0.5, cl3};
Point(11) = {0, 0, 0.5, cl3};

Circle(1) = {3, 4, 1};
Circle(2) = {1, 4, 2};
Circle(3) = {6, 4, 1};
Circle(4) = {1, 4, 5};
Circle(5) = {3, 4, 6};
Circle(6) = {6, 4, 2};
Circle(7) = {2, 4, 5};
Circle(8) = {5, 4, 3};
Line Loop(1) = {3, 2, -6};
Ruled Surface(1) = {1};
Line Loop(2) = {-2, -7, 4};
Ruled Surface(2) = {2};
Line Loop(3) = {-4, -8, -1};
Ruled Surface(3) = {3};
Line Loop(4) = {1, -3, -5};
Ruled Surface(4) = {4};
Line Loop(5) = {6, 7, 8, 5};
Plane Surface(5) = {5};
Surface Loop(1) = {3, 2, 1, 4, 5};
Volume(1) = {1};

Circle(9) = {9, 7, 11};
Circle(10) = {11, 7, 8};
Circle(11) = {8, 7, 10};
Circle(12) = {10, 7, 9};
Line Loop(6) = {10, 11, 12, 9};
Plane Surface(6) = {6};
Extrude {0, -cl3, 0} {
  Surface{6};
}

Physical Surface("top_surface") = {5};
Physical Surface("rigid_surface") = {6};
Physical Volume("top_body") = {1};
Physical Volume("bottom_body") = {2};
Physical Surface("contact_surface") = {1, 2, 4, 3};
