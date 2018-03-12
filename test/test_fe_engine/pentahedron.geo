l = .5;
Point(1) = {-.5, 0, 0, l};

Extrude {1, 0, 0} { Point{1}; }
Extrude {0, 1, 0} { Line{1}; }
Extrude {0, 0, 1} { Surface{5}; }

Transfinite Surface "*";
Transfinite Volume "*";

Recombine Surface {27, 22, 5, 14};
Mesh.SecondOrderIncomplete = 1;
