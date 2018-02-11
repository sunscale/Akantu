L = 3;
s = 1;

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

Physical Surface("left") = {26};
Physical Surface("right") = {18};
