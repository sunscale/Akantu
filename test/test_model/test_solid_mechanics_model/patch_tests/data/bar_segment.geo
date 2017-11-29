h = .6;
Point(1) = {0, 0, 0, h};
Point(2) = {10, 0, 0, h};
Line(1) = {1, 2};
Physical Point("Left") = {1};
Physical Point("Right") = {2};
//+
Physical Line("bulk") = {1};
