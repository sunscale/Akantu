Point(1) = {0, 0, 0, 1};
Point(2) = {2, 0, 0, 1};
Point(3) = {2, 1, 0, 1};
Point(4) = {0, 1, 0, 1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {4, 1, 2, 3};
Plane Surface(6) = {5};
Physical Surface(7) = {6};

Transfinite Line {1, 3} = 3;
Transfinite Line {2, 4} = 2;

Recombine Surface "*";
