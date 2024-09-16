//+
SetFactory("OpenCASCADE");
Point(1) = {-15, 0, 0, 1.0};
Point(2) = {+15, 0, 0, 1.0};
Line(1) = {1, 2};
//Rectangle(1) = {-75, -50, 0, 150, 100, 0};
Rectangle(1) = {-50, -40, 0, 100, 80, 0};
//+
BooleanFragments{ Surface{1}; Delete; }{ Curve{1}; Delete; }

MeshSize {:} = 10;
MeshSize {1,2} = 1;

// "Frontal-Delaunay" 2D meshing algorithm (Mesh.Algorithm = 6)
// usually leads to the highest quality meshes, the
// "Delaunay" algorithm (Mesh.Algorithm = 5) will handle complex mesh size


Mesh.Algorithm = 6; // Frontal-Delaunay
//Mesh.Algorithm = 9; //  Packing of Parallelograms
//Mesh.RecombinationAlgorithm = 2;
// Mesh recombination algorithm (0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad)
//Recombine Surface{:};

Mesh.SubdivisionAlgorithm = 1;
//Mesh subdivision algorithm (0: none, 1: all quadrangles, 2: all hexahedra, 3:barycentric)
//RefineMesh;


//Physical Surface("fault") = {1};


//
