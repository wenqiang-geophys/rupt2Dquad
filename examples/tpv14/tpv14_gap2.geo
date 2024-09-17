// Gmsh project created on Wed Sep 22 12:42:09 2021
//+
SetFactory("OpenCASCADE");
//+
Point(2) = {0, 0, 0, 1.0};
//+
Point(3) = {12e3, 0, 0, 1.0};
//+
Point(4) = {50e3, 0, 0, 1.0};
//+
Point(5) = {-50e3, 0, 0, 1.0};
//+
Point(6) = {-16e3, 0, 0, 1.0};
//+
Point(7) = {10.392304845413264e3, -6e3, 0, 1.0};
//+
Point(8) = { 43.301270189221931e3, -25e3, 0, 1.0};
//+
Point(9) = { 86.602540378443862, -50, 0, 1.0};
//+
Circle(1) = {5, 2, 4};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 4};

//+
Line(8) = {9, 7};
//+
Line(9) = {7, 8};
//+
Line(10) = {2, 9};
//+
Circle(11) = {4, 2, 8};
//+
Circle(12) = {8, 2, 5};
//+
Curve Loop(1) = {1, -7, -6, -5, -4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 10, 8, 9, 12, 4};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {6, 7, 11, -9, -8, -10};
//+
Plane Surface(3) = {3};
//+
MeshSize {7, 3, 6, 9, 2} = 100;
//+
MeshSize {5, 4, 8} = 1e4;
//+
MeshSize {2, 9} = 100;
//+
MeshSize {5, 4} = 1e4;
//+
MeshSize {8} = 1e4;
