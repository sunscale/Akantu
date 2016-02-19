// Concrete geometry

// Target mesh size
lc = 0.3;

Point(1) = {0, 0, 0, lc};
Point(2) = {10, 0, 0, lc};
Point(3) = {10, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Line(1) = {4, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};

Line Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

// Boundary conditions
Physical Line("XBlocked") = {3};
Physical Point("YBlocked") = {1};
Physical Line("Force") = {4};

Physical Surface("concrete") = {1};

SecondOrderIncomplete=1;