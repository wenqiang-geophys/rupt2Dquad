//+
SetFactory("OpenCASCADE");
Point(1) = {-10e3, 0, 0, 1.0};
Point(2) = {+10e3, 0, 0, 1.0};
Line(1) = {1, 2};
Lx = 40e3;
Ly = 20e3;
Rectangle(1) = {-Lx/2, -Ly/2, 0, Lx, Ly, 0};
//+
BooleanFragments{ Surface{1}; Delete; }{ Curve{1}; Delete; }

MeshSize {:} = 200;
MeshSize {1,2} = 200;

// "Frontal-Delaunay" 2D meshing algorithm (Mesh.Algorithm = 6)
// usually leads to the highest quality meshes, the
// "Delaunay" algorithm (Mesh.Algorithm = 5) will handle complex mesh size


Mesh.Algorithm = 6; // Frontal-Delaunay
//Mesh.Algorithm = 9; //  Packing of Parallelograms
Mesh.RecombinationAlgorithm = 2;
// Mesh recombination algorithm (0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad)
Recombine Surface{:};

//Mesh.SubdivisionAlgorithm = 1;
//Mesh subdivision algorithm (0: none, 1: all quadrangles, 2: all hexahedra, 3:barycentric)
//RefineMesh;


//Physical Surface("fault") = {1};


