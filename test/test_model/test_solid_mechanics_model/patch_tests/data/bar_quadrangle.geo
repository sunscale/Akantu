h = .4;
Point(1) = {0, 0, 0, h};
Point(2) = {10, 0, 0, h};

Point(3) = {0, 1., 0, h};
Point(4) = {10, 1., 0, h};

Line(1) = {1, 2};
//+
Line(2) = {4, 3};
//+
Line(3) = {2, 4};
//+
Line(4) = {3, 1};
//+
Line Loop(1) = {1, 3, 2, 4};
//+
Plane Surface(1) = {1};


Physical Line("Left") = {4};
Physical Line("Right") = {3};
//+
Physical Surface("bulk") = {1};
//+
Recombine Surface {1};
Mesh.SecondOrderIncomplete = 1;