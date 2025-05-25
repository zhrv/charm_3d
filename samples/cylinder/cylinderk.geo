// Gmsh project created on Sun May 18 00:21:41 2025
SetFactory("OpenCASCADE");
//+
Circle(1) = {-0.5, 0.5, 0, 0.1, 0, 2*Pi};
//+
Cylinder(1) = {0, 0.5, 0.5, 1, 0, 0, 0.5, 2*Pi};
//+
Rotate {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Curve{1}; 
}

////без разреза
//+
Curve Loop(4) = {4};
//+
Curve Loop(5) = {1};
//+
Plane Surface(4) = {4, 5};
//+
Plane Surface(5) = {5}; //разрез

////-------

//+
Physical Surface("LEFT_WALL", 6) = {4};
//+
Physical Surface("LEFT_INLET", 7) = {5};
//+
Physical Surface("RIGHT", 8) = {2};
//+
Physical Surface("WALL", 9) = {1};
//+
Physical Volume("CYLINDER", 10) = {1};
