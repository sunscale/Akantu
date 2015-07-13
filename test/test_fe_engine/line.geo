// Mesh size
h = 0.1;

Lx = 1;

Point(1) = { 0.0, 0.0, 0.0, h};
Point(2) = { Lx,  0.0, 0.0, h};

Line(1) = {1,2};
