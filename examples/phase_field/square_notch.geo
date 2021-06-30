element_size = 0.050;

Point(1) = {0.5, 0.5, 0, element_size};
Point(2) = {-0.5, 0.5, 0, element_size};
Point(3) = {-0.5, -0.5, 0, element_size};
Point(4) = {0.5, -0.5, 0, element_size};
Point(5) = {-0.5, 0.001, 0, element_size};
Point(6) = {0., 0.0, 0, element_size};
Point(7) = {0.5, 0.0, 0, element_size};
Point(8) = {-0.5, -0.001, 0, element_size};

Line(1) = {3, 4};
Line(2) = {4, 7};
Line(3) = {7, 1};
Line(4) = {1, 2};
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 8};
Line(8) = {8, 3};

Line Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};

Plane Surface(1) = {1};

Physical Surface("plate") = {1};

Physical Line("bottom") = {1};
Physical Line("right") = {2, 3};
Physical Line("top") = {4};
Physical Line("left") = {5,8};
