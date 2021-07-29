
N = 10;
CellSize = 0;

Point(1)={ -1000,  0, 0, CellSize};
Point(2)={ 0,  0, 0, CellSize};
Point(3)={ 1000,  176.33, 0, CellSize};
Point(4)={ 6000,  176.33, 0, CellSize};
Point(5)={ 6000,  1000, 0, CellSize};
Point(6)={ 1000,  1000, 0, CellSize};
Point(7)={ 0,  1000, 0, CellSize};
Point(8)={-1000,  1000, 0, CellSize};


Line(1)  = {1,2};
Line(2)  = {2,3};
Line(3)  = {3,4};
Line(4)  = {4,5};
Line(5)  = {6,5};
Line(6)  = {7,6};
Line(7)  = {8,7};
Line(8)  = {1,8};
Line(9)  = {2,7};
Line(10)  = {3,6};

Transfinite Line {  1, 2, 6, 7, 8, 4, 9, 10 }  = N+1   Using Progression 1;
Transfinite Line {  3, 5 }              = 6*N+1 Using Progression 1;


Line Loop(101) = {1,9,-7,-8};
Plane Surface(1) = {101};
Recombine Surface (1);
Transfinite Surface (1);

Line Loop(102) = {2,10,-6,-9};
Plane Surface(2) = {102};
Recombine Surface (2);
Transfinite Surface (2);

Line Loop(103) = {3,4,-5,-10};
Plane Surface(3) = {103};
Recombine Surface (3);
Transfinite Surface (3);


Extrude {0, 0, 100} {
  Surface{1, 2, 3, 4}; Layers{1}; Recombine;
}


Physical Volume("KLIN")          = {1,2,3};


Physical Surface("FRONT") = {125,147,169};
Physical Surface("BACK") = {1,2,3};

Physical Surface("LEFT") = {124};
Physical Surface("RIGHT") = {160};

Physical Surface("TOP") = {120,142,164};
Physical Surface("BOTTOM") = {112,134,156};
