l = 100;
Point(1) = {0, 0, 0, l};
Point(2) = {0, 1, 0, l};
Point(3) = {0, 0, 1, l};
Point(4) = {1, 0, 0, l};
Point(5) = {1, 1, 0, l};
Point(9) = {1, 0, 1, l};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line(7) = {4, 5};
Line(8) = {5, 9};
Line(9) = {9, 4};
Line(11) = {1, 4};
Line(12) = {2, 5};
Line(16) = {3, 9};
Line Loop(5) = {1, 2, 3};
Plane Surface(5) = {5};
Line Loop(13) = {1, 12, -7, -11};
Ruled Surface(13) = {13};
Line Loop(17) = {2, 16, -8, -12};
Ruled Surface(17) = {17};
Line Loop(21) = {3, 11, -9, -16};
Ruled Surface(21) = {21};
Line Loop(22) = {7, 8, 9};
Plane Surface(22) = {22};
Surface Loop(1) = {5, 22, 13, 17, 21};
Volume(1) = {1};

Rotate {{1, 0, 0}, {0, 0.5, 0.5}, Pi} {
  Duplicata { Volume{1}; }
}

Transfinite Line {16, 11, 12, 36} = 2 Using Progression 1;
Transfinite Surface "*";
Recombine Surface {17, 21, 13, 32, 42};
Transfinite Volume "*";

Mesh.SecondOrderIncomplete = 1;
