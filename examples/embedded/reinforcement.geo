// Reinforcement geometry

// Target mesh size
lc = 1;

Point(1) = {0, 0.25, 0, lc};
Point(2) = {10, 0.25, 0, lc};

Line(1) = {1, 2};

Physical Line("reinforcement") = {1};
