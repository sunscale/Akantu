h = 0.5;
length = 2;

Point(1) = {-length/2, -length/2, 0, h};
Point(2) = {0, -length/2, 0, h};
line1[] = Extrude{length/2, 0, 0}{ Point{1}; };
line2[] = Extrude{length/2, 0, 0}{ Point{2}; };
surface1[] = Extrude{0, length, 0}{ Line{line1[1]}; };
surface2[] = Extrude{0, length, 0}{ Line{line2[1]}; };

Physical Surface("outside") = {surface1[1]};
Physical Surface("inside") = {surface2[1]};

Transfinite Surface "*" AlternateLeft;
