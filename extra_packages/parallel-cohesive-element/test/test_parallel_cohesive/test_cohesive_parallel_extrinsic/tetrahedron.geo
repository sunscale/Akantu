lc = 1;
h = 1;

Point(1) = {lc, lc, lc, h};

line[] = Extrude{-2*lc,0,0}{ Point{1}; };
surface[] = Extrude{0,-2*lc,0}{ Line{line[1]}; };
volume[] = Extrude{0,0,-2*lc}{ Surface{surface[1]}; };

Physical Volume(1) = {volume[1]};

Transfinite Surface "*";
Transfinite Volume "*";
