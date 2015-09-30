L = 1;

Point(1) = {-L/2., 0, 0, L/2};

line[] = Extrude{L,0,0}{ Point{1}; };
surface[] = Extrude{0,L,0}{ Line{line[1]}; };
volume[] = Extrude{0,0,L/2}{ Surface{surface[1]}; };

Physical Volume(1) = {volume[1]};
//Physical Surface(1) = {surface[1]};

Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
