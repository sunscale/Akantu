h = 0.7;

z = 1;
y = 2;
x = 2;

Point(1) = { x/2, y/2, -z/2, h};
Point(2) = {-x/2, y/2, -z/2, h};
Point(3) = {-x/2,-y/2, -z/2, h};
Point(4) = { x/2,-y/2, -z/2, h};

Point(5) = {-x/2, 0,  -z/2, h};
Point(6) = { x/2, 0,  -z/2, h};

Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 6};
Line(4) = {6, 1};

Line(5) = {5, 3};
Line(6) = {3, 4};
Line(7) = {4, 6};

Line Loop(1) = {1, 2, 3, 4};
Line Loop(2) = {-3, 5, 6, 7};
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Extrude {0, 0, z} {
  Surface{1}; Surface{2};
}

Physical Surface("fixed") = {46};
Physical Surface("insertion") = {24};
Physical Surface("loading") = {16};
Physical Surface("sides") = {1, 20, 29, 28, 50, 2, 51, 42};

Physical Volume("body") = {1, 2};
