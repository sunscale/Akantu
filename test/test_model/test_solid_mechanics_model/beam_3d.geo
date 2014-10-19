cl1 = 0.8;
Point(1)  = {0, 0, 0, cl1};
Point(2)  = {0, 1, 0, cl1};
Point(3)  = {0, 1, 1, cl1};
Point(4)  = {0, 0, 1, cl1};
Point(5)  = {0, 0.5, 1, cl1};
Point(6)  = {0, 0.5, 0, cl1};
Point(7)  = {0, 0.5, 0.5, cl1};

Line(1) = {1, 4};
Line(2) = {4, 5};
Line(3) = {5, 3};
Line(4) = {3, 2};
Line(5) = {2, 6};
Line(6) = {6, 1};
Line(7) = {6, 7};
Line(8) = {7, 5};
Line Loop(9) = {6, 1, 2, -8, -7};
Plane Surface(10) = {9};
Line Loop(11) = {3, 4, 5, 7, 8};
Plane Surface(12) = {11};
Extrude {5, 0, 0} {
  Surface{12, 10};
}
Extrude {5, 0, 0} {
  Surface{39, 66};
}

Delete {
  Volume{3, 4, 1, 2};
}
Delete {
  Surface{38, 34, 88, 92};
}
Delete {
  Line{87, 83, 74, 20, 33, 29};
}

Surface Loop(121) = {12, 22, 26, 30, 10, 49, 53, 57, 39, 66};
Volume(122) = {121};
Surface Loop(123) = {103, 107, 111, 120, 66, 39, 76, 80, 84, 93};
Volume(124) = {123};
Delete {
  Volume{124, 122};
}
Delete {
  Surface{84, 103, 76, 111, 57, 22, 30, 49};
}
Delete {
  Line{83, 74, 20, 29};
}
Line Loop(124) = {79, 70, 95, -102, -41, -16};
Plane Surface(125) = {124};
Line Loop(126) = {16, 41, -48, -6, -5, 25};
Plane Surface(127) = {126};
Line Loop(128) = {21, -14, -43, -52, 2, 3};
Plane Surface(129) = {128};
Line Loop(130) = {75, -68, -97, -106, 43, 14};
Plane Surface(131) = {130};
Surface Loop(132) = {131, 80, 125, 93, 120, 107, 39, 66};
Volume(133) = {132};
Surface Loop(134) = {26, 12, 129, 53, 10, 127, 39, 66};
Volume(135) = {134};
Physical Volume(136) = {133, 135};
