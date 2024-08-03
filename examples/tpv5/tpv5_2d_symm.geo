//+
SetFactory("OpenCASCADE");

faultsize = 0.5;
size = 20;

Point(1) = {-15, 0, 0, 1.0};
Point(2) = {+15, 0, 0, 1.0};
Line(1) = {1, 2};
Point(3) = {-15, -2*faultsize, 0, 1.0};
Point(4) = {+15, -2*faultsize, 0, 1.0};
Line(2) = {3, 4};
Point(5) = {-15, +2*faultsize, 0, 1.0};
Point(6) = {+15, +2*faultsize, 0, 1.0};
Line(3) = {5, 6};
Point(7) = {-15, -4*faultsize, 0, 1.0};
Point(8) = {+15, -4*faultsize, 0, 1.0};
Line(4) = {7, 8};
Point(9) = {-15, +4*faultsize, 0, 1.0};
Point(10)= {+15, +4*faultsize, 0, 1.0};
Line(5) = {9, 10};
//Rectangle(1) = {-75, -50, 0, 150, 100, 0};
Rectangle(1) = {-50, -40, 0, 100, 80, 0};
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
