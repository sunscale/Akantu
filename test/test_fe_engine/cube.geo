
Point(1) = {-0.5, -0.5, -0.5, 1.0};
l[] = Extrude {1, 0, 0} {
  Point{1}; Layers {2};
};
s[] = Extrude {0, 1, 0} {
  Line{l[1]}; Layers {2}; Recombine;
};
Extrude {0, 0, 1} {
  Surface{s[1]}; Layers{2}; Recombine;
}


Mesh.SecondOrderIncomplete = 1;
