// Gmsh project created on Wed Sep 22 12:42:09 2021

SetFactory("OpenCASCADE");
faultsize = 200*2;
boundsize = 5e3*2;
Point(1) = {0, 0, 0, 1.0};
Point(2) = {12e3, 0, 0, 1.0};
Point(3) = {-16e3, 0, 0, 1.0};

theta = -30;
theta = theta/180*Pi;
Point(4) = {12e3*Cos(theta), 12e3*Sin(theta), 0, 1.0};
Line(1) = {3, 1};
Line(2) = {1, 2};
Line(3) = {1, 4};
Circle(4) = {0, 0, 0, 50e3, 0, 2*Pi};
Curve Loop(1) = {4};
Surface(1) = {1};
BooleanFragments{ Surface{1}; Delete; }{ Curve{1:3}; Delete; }
//+
MeshSize {:} = boundsize;
MeshSize {1:4} = faultsize;
//+
//Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.Algorithm = 9; //  Packing of Parallelograms
Mesh.RecombinationAlgorithm = 2;
Mesh.RecombineAll = 1;

Mesh.RecombineAll = 0;
Mesh.Algorithm = 5; // Delaunay
Mesh.SubdivisionAlgorithm = 1;
