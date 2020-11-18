h = .5;
length = 2;

Point(1) = {-length/2, -length/2, 0, h};
Point(2) = {0, -length/2, 0, h};
line1[] = Extrude{length/2, 0, 0}{ Point{1}; };
line2[] = Extrude{length/2, 0, 0}{ Point{2}; };
surface1[] = Extrude{0, length, 0}{ Line{line1[1]}; };
surface2[] = Extrude{0, length, 0}{ Line{line2[1]}; };

Transfinite Surface "*";
Recombine Surface "*";

Physical Surface("mat_1") = {surface1[1]};
Physical Surface("mat_2") = {surface2[1]};
