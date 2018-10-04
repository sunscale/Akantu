element_size = 1;

Point(1) = {0.5, 0.5, 0, element_size};
Point(2) = {-0.5, 0.5, 0, element_size};
Point(3) = {-0.5, -0.5, 0, element_size};
Point(4) = {0.5, -0.5, 0, element_size};

Line(1) = {3, 4};
Line(2) = {4, 1};
Line(3) = {1, 2};
Line(4) = {2, 3};

Line Loop(1) = {1, 2, 3, 4};

Plane Surface(1) = {1};

Physical Line("bottom") = {1};
Physical Line("right") = {2};
Physical Line("top") = {3};
Physical Line("left") = {4};

Physical Surface("plate") = {1};
