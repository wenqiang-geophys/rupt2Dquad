//+
SetFactory("OpenCASCADE");
Rectangle(1) = {-50, 0, 0, 100, -100, 0};
size = 2;
MeshSize {:} = size;

Mesh.Algorithm = 8; // Frontal-Delaunay for quads
//Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.RecombinationAlgorithm = 2; // 2 or 3
Recombine Surface{:};
