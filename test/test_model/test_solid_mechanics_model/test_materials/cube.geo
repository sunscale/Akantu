Point(1) = {-.5, -.5, -.5, 0.1};

Extrude {1, 0, 0} {
  Point{1}; 
}

Extrude {0, 1, 0} {
  Line{1}; 
}

Line Loop(1) = {3, 2, -4, -1};
Plane Surface(6) = {1};

Extrude {0, 0, 1} {
  Surface{5}; 
}

Physical Surface("box") = {15, 27, 23, 28, 19, 5};
Physical Volume("bulk") = {1};
