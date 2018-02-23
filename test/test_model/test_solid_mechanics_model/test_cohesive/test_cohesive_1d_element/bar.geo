l = 1;
n = 4;

Point(1) = {-l/2., 0, 0, 1.0};
Point(2) = {l/2., 0, 0, 1.0};
Line(1) = {1, 2};
Physical Line(2) = {1};
Physical Point ("left") = {1};
Physical Point ("right") = {1};
Transfinite Line {1} = n + 1;
