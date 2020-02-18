h = 1.;

Point(1) = { 1,  1, 0, h};
Point(2) = {-1,  1, 0, h};
Point(3) = {-1, -1, 0, h};
Point(4) = { 1, -1, 0, h};
Point(5) = {-1,  0, 0, h};
Point(6) = { 1,  0, 0, h};

Line(1) = {3, 4};
Line(2) = {4, 6};
Line(3) = {6, 5};
Line(4) = {5, 3};
Line(5) = {6, 1};
Line(6) = {1, 2};
Line(7) = {2, 5};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Curve Loop(2) = {-3, 5, 6, 7};
Plane Surface(2) = {2};

Physical Curve("top") = {6};
Physical Curve("bottom") = {1};

Physical Curve("left") = {4, 7};
Physical Curve("right") = {2, 5};

Physical Surface("btop") = {2};
Physical Surface("bbottom") = {1};
