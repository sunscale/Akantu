L = 3;
s = 3;

Point(1) = {-L/2., -L/2., 0, s};

line[] = Extrude{L, 0, 0} {
  Point{1};
};

surface[] = Extrude{0, L, 0} {
  Line{line[1]};
};

volume[] = Extrude{0, 0, L} {
  Surface{surface[1]};
};

Physical Volume("body") = {volume[1]};

Periodic Surface 26 {12, -10, -21, -3} = 18 {13, 8, -17, -4};
Periodic Surface 5 {1, 4, -2, -3} = 27 {7, 8, 9, 10};
Periodic Surface 14 {13, 1, -12, -7} = 22 {17, 2,  -21, 9};

Physical Surface("left") = {26};
Physical Surface("right") = {18};

Physical Surface("front") = {27};
Physical Surface("back") = {-5};

Physical Surface("bottom") = {22};
Physical Surface("top") = {14};
