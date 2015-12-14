dx = 0.01;
w = 0.2;
radius = 0.015;

Point(1) = {0,0,0,dx};
Point(2) = {0,.03,0,dx};
Point(3) = {0,-.03,0,dx};
Point(4) = {w,0,0,dx};
Point(5) = {w,.03,0,dx};
Point(6) = {w,-.03,0,dx};
Point(7) = {0,w,0,dx};
Point(8) = {0,w,0.3,dx};
Point(9) = {0,w,-.03,dx};
Point(10) = {w,w,0,dx};
Point(11) = {w,w,.03,dx};
Point(12) = {w,w,-.03,dx};
Point(13) = {0.5*w,0.5*w,0,dx};
Point(14) = {0.5*w+radius,0.5*w+radius,0,dx};
Point(15) = {0.5*w-radius,0.5*w+radius,0,dx};
Point(16) = {0.5*w+radius,0.5*w-radius,0,dx};
Point(17) = {0.5*w-radius,0.5*w-radius,0,dx};
Line(1) = {1, 2};
Line(2) = {2, 5};
Line(3) = {5, 4};
Line(4) = {1, 4};
Line(5) = {1, 3};
Line(6) = {6, 4};
Line(7) = {3, 6};
Line(8) = {8, 11};
Line(9) = {7, 10};
Line(10) = {9, 12};
Line(11) = {8, 7};
Line(12) = {11, 10};
Line(13) = {10, 12};
Line(14) = {7, 9};
Line(15) = {2, 8};
Line(16) = {1, 7};
Line(17) = {3, 9};
Line(18) = {5, 11};
Line(19) = {4, 10};
Line(20) = {6, 12};
Circle(41) = {16, 13, 14};
Circle(42) = {14, 13, 15};
Circle(43) = {15, 13, 17};
Circle(44) = {17, 13, 16};
Translate {0.05, 0.05, 0} {
  Duplicata { Line{44, 43, 42, 41}; }
}
Translate {-0.05, 0.05, 0} {
  Duplicata { Line{44, 43, 42, 41}; }
}
Translate {0.05, -0.05,0} {
  Duplicata { Line{44, 43, 42, 41}; }
}
Translate {-0.05,-0.05,0} {
  Duplicata { Line{44, 43, 42, 41}; }
}
Delete {
  Line{2, 15, 1, 3, 18, 8, 7, 5, 17, 14, 10, 13, 12, 20};
}
Delete {
  Line{11, 6};
}
Delete {
  Point{2, 3, 9, 12, 11, 5, 8, 6};
}
Line Loop(61) = {4, 19, -9, -16};
Line Loop(62) = {60, 59, 58, 57};
Line Loop(63) = {44, 41, 42, 43};
Line Loop(64) = {55, 54, 53, 56};
Line Loop(65) = {52, 51, 50, 49};
Line Loop(66) = {46, 45, 48, 47};
Plane Surface(67) = {61, 62, 63, 64, 65, 66};
Physical Surface("cheese") = {67};
Physical Line("holes") = {60, 57, 58, 59, 56, 53, 54, 55, 41, 44, 43, 42, 52, 51, 50, 49, 45, 48, 47, 46};
Physical Line("crust") = {4, 19, 9, 16};
