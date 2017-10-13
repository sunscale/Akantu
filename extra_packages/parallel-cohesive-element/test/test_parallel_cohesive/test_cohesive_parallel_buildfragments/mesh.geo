L = 0.03;
h = L/5;
s = L/20;

Point(1) = {-L/2., 0, 0, s};

line[] = Extrude{L,0,0}{ Point{1}; };
surface[] = Extrude{0,h,0}{ Line{line[1]}; };
volume[] = Extrude{0,0,h}{ Surface{surface[1]}; };

Physical Volume(1) = {volume[1]};

Transfinite Surface "*";
Transfinite Volume "*";
