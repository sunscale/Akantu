
L = 10;
l = 0.1;
h = l;

Point(1) = {0, 0, 0, h};
Point(2) = {L, 0, 0, h};
Point(3) = {L, l, 0, h};
Point(4) = {0, l, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Surface "*";


Physical Surface("bulk") = {1};
Physical Line("left") = {4};
Physical Line("bottom") = {1};
Physical Line("top") = {3};
Physical Line("right") = {2};
