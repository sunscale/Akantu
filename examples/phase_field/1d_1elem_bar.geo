lc = 1./16;
Point(1) = {-0.5, 0, 0, lc};
Point(2) = {0.5, 0, 0, lc};
Line(1) = {1, 2};
Physical Line(2) = {1};

//Physical names
Physical Point("Left") = {1};
Physical Point("Right") = {2};

