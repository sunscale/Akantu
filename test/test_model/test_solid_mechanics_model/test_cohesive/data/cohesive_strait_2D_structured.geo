h = 0.2;

Point(1) = { 1, 1, 0, h};
Point(2) = {-1, 1, 0, h};
Point(3) = {-1,-1, 0, h};
Point(4) = { 1,-1, 0, h};

Point(5) = {-1, 0, 0, h};
Point(6) = { 1, 0, 0, h};

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

Physical Line("fixed") = {6};
Physical Line("loading") = {1};
Physical Line("insertion") = {3};
Physical Line("sides") = {2, 5, 7, 4};

Physical Surface("body") = {1, 2};

Recombine Surface "*";
Transfinite Surface "*";
Mesh.SecondOrderIncomplete = 1;
