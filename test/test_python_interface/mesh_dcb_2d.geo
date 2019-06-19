//Mesh size
dx = 57.8e-5;

Point(1) = {0,0,0,dx};
Point(2) = {0,0.00055,0,dx};
Point(3) = {0,-0.00055,0,dx};
Point(4) = {57.8e-3,0,0,dx};
Point(5) = {57.8e-3,0.00055,0,dx};
Point(6) = {57.8e-3,-0.00055,0,dx};
Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 4};
Line(4) = {1, 4};
Line(5) = {1, 3};
Line(6) = {6, 4};
Line(7) = {3, 6};
Line Loop(8) = {2, 3, -4, 1};
Plane Surface(9) = {-8};
Line Loop(10) = {5, 7, 6, -4};
Plane Surface(11) = {10};
Physical Surface("bulk") = {9,11};
Physical Line("coh") = {4};
Physical Line("edge") = {1};
Transfinite Surface "*";
Recombine Surface "*";
Mesh.SecondOrderIncomplete = 1;
