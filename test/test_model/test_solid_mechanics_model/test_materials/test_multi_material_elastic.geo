// Mesh size
h  = 0.5;

// Dimensions of the bar
Lx = .5;
Ly = .5;

// Geometry
Point(1) = { 0.0, 0.0, 0.0, h};
Point(2) = { Lx,  0.0, 0.0, h};

Line(3) = {1, 2};

Extrude {0, Ly, 0} {
  Line{3};
}

Extrude {0, Ly, 0} {
  Line{4};
}

Extrude {0, Ly, 0} {
  Line{8};
}

Physical Point("corner") = {1};

Physical Line("ground") = {3};
Physical Line("air") = {12};

Physical Surface("strong") = {11};
Physical Surface("weak") = {7, 15};
