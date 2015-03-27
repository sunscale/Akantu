f = 1.0;
Point(1) = {0, 0, 0, f};
Point(2) = {1, 0, 0, f};
Point(3) = {1, 1, 0, f};
Point(4) = {0, 1, 0, f};
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {1, 2};
Line(4) = {4, 1};
Line Loop(5) = {4, 3, -2, -1};
Plane Surface(6) = {5};
Extrude {0, 0, 1} {
  Surface{6};
}
