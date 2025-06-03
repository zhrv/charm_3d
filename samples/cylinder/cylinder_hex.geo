Mesh.MshFileVersion=2.2;
SetFactory("OpenCascade");

// Центр круга
Point(1) = {0, 0.5, 0.5, 1.0};

// Внешние точки круга (радиус 0.5)
Point(2) = {0.5, 0.5, 0.5, 1.0};   // Правая
Point(3) = {0, 1.0, 0.5, 1.0};     // Верхняя
Point(4) = {-0.5, 0.5, 0.5, 1.0};  // Левая
Point(5) = {0, 0.0, 0.5, 1.0};     // Нижняя

// Внутренние точки (радиус 0.1)
Point(6) = {0.1, 0.5, 0.5, 1.0};   // Правая
Point(7) = {0, 0.6, 0.5, 1.0};     // Верхняя
Point(8) = {-0.1, 0.5, 0.5, 1.0};  // Левая
Point(9) = {0, 0.4, 0.5, 1.0};     // Нижняя

Translate {-0.5, 0, 0} {
  Point{1:9};
}

Translate {0, 0, -0.5} {
  Point{1:9};
}


// Внешний круг
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Внутреннее отверстие
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 6};
Rotate {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Curve{1}; Curve{2}; Curve{3}; Curve{4}; Curve{5}; Curve{6}; Curve{7}; Curve{8}; 
}

// Соединительные линии для деления на 4 части
Line(9) = {2, 6};  // Правая
Line(10) = {3, 7}; // Верхняя
Line(11) = {4, 8}; // Левая
Line(12) = {5, 9}; // Нижняя

// 4 замкнутых сектора
Line Loop(21) = {1, 10, -5, -9};
Plane Surface(31) = {21};

Line Loop(22) = {2, 11, -6, -10};
Plane Surface(32) = {22};

Line Loop(23) = {3, 12, -7, -11};
Plane Surface(33) = {23};

Line Loop(24) = {4, 9, -8, -12};
Plane Surface(34) = {24};
/////////////////////////////////
Recombine Surface {31};
Recombine Surface {32};
Recombine Surface {33};
Recombine Surface {34};
//+
Transfinite Surface {34};
//+
Transfinite Surface {31};
//+
Transfinite Surface {32};
//+
Transfinite Surface {33};
//+
Transfinite Curve {8, 4} = 15 Using Progression 1;
//+
Transfinite Curve {9, 12} = 15 Using Progression 1;
//+
Transfinite Curve {5, 1} = 15 Using Progression 1;
//+
Transfinite Curve {10, 11} = 15 Using Progression 1;
//+
Transfinite Curve {6, 2} = 15 Using Progression 1;
//+
Transfinite Curve {7, 3} = 15 Using Progression 1;

//отдельно то же самое для разреза

//+
Curve Loop(25) = {8, 5, 6, 7};
//+
Plane Surface(35) = {25};
//+
Transfinite Surface {35};
//+
Transfinite Curve {5, 6, 7, 8} = 15 Using Progression 1;
//+
Recombine Surface {35};



//+
Extrude {1, 0, 0} {
  Surface{32}; Surface{33}; Surface{34}; Surface{31}; Surface{35}; Layers {10}; Recombine;
}


//+
Physical Surface("LEFT_WALL", 146) = {34, 33, 32, 31};
//+
Physical Surface("LEFT_INLET", 147) = {35};
//+
Physical Surface("RIGHT", 148) = {57, 79, 101, 123, 145};
//+
Physical Surface("WALL", 149) = {110, 44, 66, 88, 118};
//+
Physical Volume("CYLINDER", 150) = {4, 1, 5, 3, 2};
