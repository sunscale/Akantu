//+
Point(1) = {-0.9, -0.5, -0, 1.0};
//+
Point(2) = {0.3, -0.5, 0, 1.0};
//+
Point(3) = {0.6, 0.2, 0, 1.0};
//+
Point(4) = {0.2, 0.8, 0, 1.0};
//+
Point(5) = {-0.7, 0.7, -0, 1.0};
//+
Point(6) = {-1.1, 0, -0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Line(7) = {2, 5};
//+
Line Loop(8) = {1, 7, 5, 6};
//+
Plane Surface(9) = {8};
//+
Line Loop(10) = {7, -4, -3, -2};
//+
Plane Surface(11) = {10};
//+
Extrude {0, 0, 1} {
  Surface{9};
}
//+
Extrude {0, 0, 1} {
  Surface{11};
}
//+
Physical Point("points") = {6, 1};
//+
Physical Line("lines") = {7, 23, 14, 19};
//+
Physical Surface("initial_surfaces") = {9, 11};
//+
Physical Surface("extruded_surfaces") = {33, 55};
//+
Physical Volume("v1") = {1};
//+
Physical Volume("v2") = {2};
