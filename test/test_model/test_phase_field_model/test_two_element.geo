h =  2;

Point(1) = {-1, -1, 0, h};
Point(2) = {0, -1, 0, h};
Point(3) = {0, 0, 0, h};
Point(4) = {-1, 0, 0, h};

Point(7) = {0, 1, 0, h};
Point(8) = {-1, 1, 0, h};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(7) = {3, 7};
Line(8) = {7, 8};
Line(9) = {8, 4};

Line Loop(1) = {3, 4, 1, 2};
Line Loop(2) = {7, 8, 9, -3};


Plane Surface(1) = {1};
Physical Surface("plate") = {1};


Plane Surface(2) = {2};
Physical Surface("rigid") = {2};

Transfinite Line {1, 2, 3, 4, 7, 8, 9} = 2 Using Progression 1;

//Transfinite Surface {1} = {1, 2, 3, 4};  
//Transfinite Surface {2} = {7, 8, 9, -3};

Recombine Surface{1, 2};g