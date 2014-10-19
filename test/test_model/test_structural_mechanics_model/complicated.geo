lz = 1;

Point(1) = {0, 0,0, lz};
Point(2) = {10, 0,0, lz};
Point(3) = {20, 0,0, lz};
Point(4) = {30, 0,0, lz};
Point(5) = {5, -8.66,0, lz};
Point(6) = {25, -8.66,0, lz};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {3, 6};
Line(5) = {5, 2};
Physical Line(1) = {1, 2, 3};
Physical Line(2) = {4, 5};

