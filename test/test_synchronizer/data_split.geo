a = DefineNumber[ .5, Name "Parameters/a" ];
b = DefineNumber[ .1, Name "Parameters/b" ];
lc = DefineNumber[ b, Name "Parameters/lc" ];

Point(1) = {-a, -a, -a, lc};
Extrude {2*a, 0, 0} { Point{1}; }
Extrude {0, 0, 2*a} { Line{1}; }
Extrude {0, 2*a, 0} { Surface{5}; Recombine; }

Point(15) = {-b, -b, -b, lc};
Extrude {2*b, 0, 0} { Point{15}; }
Extrude {0, 0, 2*b} { Line{28}; }
Extrude {0, 2*b, 0} { Surface{32}; Recombine; }

Surface Loop(55) = {22, 5, 14, 18, 27, 26};
Surface Loop(56) = {49, 32, 41, 45, 54, 53};

Volume(57) = {55, 56};
Physical Volume("volume") = {57};

Physical Surface("inside") = {32, 53, 41, 54, 45, 49};
Physical Surface("outside") = {27, 14, 26, 22, 18, 5};

#Transfinite Surface "*";
#Recombine   Surface "*";

#Transfinite Volume  "*";


