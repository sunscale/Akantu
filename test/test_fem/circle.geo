h  = 0.1;

Point(1) = {1.0, 0.0, 0.0,h};
Point(2) = {0.0, 1.0, 0.0,h};
Point(3) = {-1.0, 0.0, 0.0,h};
Point(4) = {0.0, -1.0, 0.0,h};
Point(5) = {0.0, 0.0, 0.0,h};

Circle(1) = {1, 5, 2};
Circle(2) = {2, 5, 3};
Circle(3) = {3, 5, 4};
Circle(4) = {4, 5, 1};

Line Loop(5) = {2, 3, 4, 1};

Plane Surface(6) = {5};
