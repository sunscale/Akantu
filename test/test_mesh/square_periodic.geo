L = 3;
s = 0.5;

Point(1) = {-L/2., -L/2., 0, s};

line[] = Extrude{L, 0, 0} {
  Point{1};
};

surface[] = Extrude{0, L, 0} {
  Line{line[1]};
};

Periodic Line {1} = {2};
Periodic Line {3} = {4};
