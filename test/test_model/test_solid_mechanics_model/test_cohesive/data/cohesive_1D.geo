h = 0.5;

Point(1) = { 1., 0., 0., h};
Point(3) = { 0., 0., 0., h};
Point(2) = {-1., 0., 0., h};

Line(1) = {2, 3};
Line(2) = {3, 1};

Physical Point ("loading") = {1};
Physical Point ("fixed") = {2};

Physical Point ("insertion") = {3};

Physical Line ("body") = {1, 2};
