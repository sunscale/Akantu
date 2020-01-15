lz = 0.3;
Point (1) = {0.00, 0, 0, lz};
Point (2) = {10, 0, 0, lz};
Point (3) = {10, 4, 0,lz};
Point (4) = {0.0 , 4,0, lz};

Line (1) = {1, 2};
Line (2) = {2, 3};
Line (3) = {3, 4};
Line (4) = {4, 1};
Line Loop (6) = {1,2, 3, 4};



Point(5) = {5, 2, 0,lz};
Point(6) = {3.5, 2, 0,lz};
Point(8) = {6.5, 2, 0,lz};
Point(9) = {5, 0.5, 0,lz};
Point(11) = {5, 3.5, 0,lz};
Point(12) = {5, 3.5, 0,lz};
Circle(7) = {6, 5, 11};
Circle(8) = {11, 5, 8};
Circle(9) = {8, 5, 9};
Circle(10) = {9, 5, 6};
Line Loop(11) = {10, 9, 8, 7};
Plane Surface(11)={6,11};
Physical Surface("Interior") = {11};
Physical Line("Traction") = {2};
Physical Line("Fixed_x") = {4};
Physical Line("Fixed_y") = {1, 3};
Physical Line("Free") = {7, 8, 9, 10};
