h = .6;
Point(1) = {0, 0, 0, h};
Point(2) = {10, 0, 0, h};

Point(3) = {0, 1., 0, h};
Point(4) = {10, 1., 0, h};

Point(5) = {0, 0, 1., h};
Point(6) = {10, 0, 1., h};

Point(7) = {0, 1., 1., h};
Point(8) = {10, 1., 1., h};



Mesh.SecondOrderIncomplete = 1;//+
Line(1) = {7, 8};
//+
Line(2) = {8, 6};
//+
Line(3) = {6, 5};
//+
Line(4) = {5, 7};
//+
Line(5) = {3, 7};
//+
Line(6) = {3, 1};
//+
Line(7) = {1, 5};
//+
Line(8) = {4, 2};
//+
Line(9) = {4, 8};
//+
Line(10) = {2, 6};
//+
Line(11) = {4, 3};
//+
Line(12) = {1, 2};
//+
Line Loop(1) = {2, -10, -8, 9};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {7, 4, -5, 6};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {12, 10, 3, -7};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {3, 4, 1, 2};
//+
Plane Surface(4) = {4};
//+
Line Loop(5) = {1, -9, 11, 5};
//+
Plane Surface(5) = {5};
//+
Line Loop(6) = {11, 6, 12, -8};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {3, 6, 5, 4, 2, 1};
//+
Volume(1) = {1};
//+
Physical Surface("BC") = {1, 2, 3, 4, 5, 6};
Physical Volume("bulk") = {1};
