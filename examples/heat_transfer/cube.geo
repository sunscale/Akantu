h = 0.1;

Point(1) = {0.0, 0.0, 0.0, h};
Point(2) = {1.0, 0.0, 0.0, h};
Point(3) = {0.0, 1.0, 0.0, h};
Point(4) = {1.0, 1.0, 0.0, h};
Point(5) = {0.0, 0.0, 1.0, h};
Point(6) = {1.0, 0.0, 1.0, h};
Point(7) = {0.0, 1.0, 1.0, h};
Point(8) = {1.0, 1.0, 1.0, h};


Line(1) = {7, 8};
Line(2) = {8, 4};
Line(3) = {4, 3};
Line(4) = {3, 7};
Line(5) = {1, 5};
Line(6) = {5, 6};
Line(7) = {6, 2};
Line(8) = {2, 1};
Line(9) = {3, 1};
Line(10) = {7, 5};
Line(11) = {8, 6};
Line(12) = {4, 2};


Line Loop(13) = {1, 11, -6, -10};
Plane Surface(14) = {13};

Line Loop(15) = {3, 4, 1, 2};
Plane Surface(16) = {15};

Line Loop(17) = {6, 7, 8, 5};
Plane Surface(18) = {17};

Line Loop(19) = {3, 9, -8, -12};
Plane Surface(20) = {19};

Line Loop(21) = {4, 10, -5, -9};
Plane Surface(22) = {21};

Line Loop(23) = {11, 7, -12, -2};
Plane Surface(24) = {23};


Surface Loop(25) = {24, 14, 16, 20, 22, 18};


Volume(26) = {25};


