el_size = 10.0;
box_size = 10.0; // this is the edge length of the box

Point(1) = {0, 0, 0, el_size};
Point(2) = {-box_size, 0, 0, el_size};
Point(3) = {-box_size, -box_size, 0, el_size};
Point(4) = {0, -box_size, 0, el_size};
Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Transfinite Surface "*";

Physical Surface(7) = {6};
Physical Line("bottom") = {4};
Physical Line("right") = {1};
Physical Line("top") = {2};
Physical Line("left") = {3};
