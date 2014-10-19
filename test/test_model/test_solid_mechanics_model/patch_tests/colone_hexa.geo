l = 10;
h = 40;
k = 2;

Point(1) = {0, 0, 0, l};

Extrude {10,0,0} {
  Point{1}; Layers{k};
}

Extrude {0,10,0} {
  Line{1}; Layers{k}; Recombine;
}

Extrude {0,0,h} {
  Surface{5}; Layers{h/(l/k)}; Recombine;
}