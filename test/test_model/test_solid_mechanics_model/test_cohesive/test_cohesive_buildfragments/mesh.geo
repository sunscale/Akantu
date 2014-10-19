L = 0.03;
s = 0.03/20.;
h = s;

Point(1) = {-L/2., 0, 0, s};
Point(2) = {L/2., 0, 0, s};
Point(3) = {L/2., h, 0, s};
Point(4) = {-L/2., h, 0, s};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Physical Surface(7) = {6};
Physical Line("Left_side") = {4};
Physical Line("Right_side") = {2};

Transfinite Line {2, 4} = 2;
Transfinite Surface "*";
Recombine Surface "*";
Mesh.SecondOrderIncomplete = 1;
