//+
SetFactory("OpenCASCADE");
Rectangle(1) = {-100, 0, 0, 200, -100, 0};
seabed = -6;
Point(5) = {-100, seabed, 0, 1.0};
Point(6) = {+100, seabed, 0, 1.0};
Line(5) = {5, 6};
BooleanFragments{ Surface{1}; Delete; }{ Curve{5}; Delete; }

size=1.5;
MeshSize {:} = size;
MeshSize {5, 6} = size*5;

Mesh.Algorithm = 8; // Frontal-Delaunay for quads
//Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.RecombinationAlgorithm = 2; // 2 or 3
Recombine Surface{:};
