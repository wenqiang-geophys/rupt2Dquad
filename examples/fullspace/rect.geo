//+
SetFactory("OpenCASCADE");
Rectangle(1) = {-100, 0, 0, 200, -50, 0};
size = 2;
MeshSize {:} = size;

Mesh.Algorithm = 8; // Frontal-Delaunay for quads
//Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.RecombinationAlgorithm = 2; // 2 or 3
Recombine Surface{:};
