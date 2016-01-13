h  = 0.1;

Lx = 10.0;
Ly = 1.0;


Point(1) = { 0.0, 0.0, 0.0, h};
Point(2) = { Lx,  0.0, 0.0, h};
Point(3) = { Lx,  Ly, 0.0,  h};
Point(4) = { 0.0, Ly, 0.0,  h};

// Base Cube
Line(12) = {1,2};
Line(23) = {2,3};
Line(34) = {3,4};
Line(41) = {4,1};

// Base Cube
Line Loop(1234) = {12,23,34,41};
Plane Surface(0) = {1234};

Physical Surface("bulk") = {0};

//Physical Point ("left1") = {4};
//Physical Point ("left2") = {1};
Physical Line ("left") = {41};
Physical Line ("right") = {23};
Physical Line ("bottom") = {12};

//Physical Point ("right1") = {2};
//Physical Point ("right2") = {3};

//Physical Point ("bottom1") = {1};
//Physical Point ("bottom2") = {2};
