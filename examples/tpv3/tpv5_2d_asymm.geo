//+
SetFactory("OpenCASCADE");

faultsize = 200;
size = 2e3;

L = 15e3;
Point(1) = {-L, 0, 0, 1.0};
Point(2) = {+L, 0, 0, 1.0};
Line(1) = {1, 2};
Point(3) = {-L, -2*faultsize, 0, 1.0};
Point(4) = {+L, -2*faultsize, 0, 1.0};
Line(2) = {3, 4};
Point(5) = {-L, +2.5*faultsize, 0, 1.0};
Point(6) = {+L, +2.5*faultsize, 0, 1.0};
Line(3) = {5, 6};
Point(7) = {-L, -4*faultsize, 0, 1.0};
Point(8) = {+L, -4*faultsize, 0, 1.0};
Line(4) = {7, 8};
Point(9) = {-L, +5*faultsize, 0, 1.0};
Point(10)= {+L, +5*faultsize, 0, 1.0};
Line(5) = {9, 10};
Rectangle(1) = {-30e3, -15e3, 0, 60e3, 30e3, 0};
//+
BooleanFragments{ Surface{1}; Delete; }{ Curve{1,2,3,4,5}; Delete; }

MeshSize {:} = size;
MeshSize {1:10} = faultsize;

// "Frontal-Delaunay" 2D meshing algorithm (Mesh.Algorithm = 6)
// usually leads to the highest quality meshes, the
// "Delaunay" algorithm (Mesh.Algorithm = 5) will handle complex mesh size


Mesh.Algorithm = 6; // Frontal-Delaunay for 2D meshes
Mesh.RecombinationAlgorithm = 2; // 2 or 3
Recombine Surface{:};
