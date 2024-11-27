//+
SetFactory("OpenCASCADE");
Point(1) = {-40e3, 0, 0, 1.0};
Point(2) = {+40e3, 0, 0, 1.0};
Line(1) = {1, 2};
Lx = 160e3;
Ly = 120e3;
Rectangle(1) = {-Lx/2, -Ly/2, 0, Lx, Ly, 0};
//+
BooleanFragments{ Surface{1}; Delete; }{ Curve{1}; Delete; }

faultsize = 200*2;
boundsize = 10e3*2;

MeshSize {:} = boundsize;
MeshSize {1,2} = faultsize;

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


Field[1] = Distance;
Field[1].CurvesList = {1};

// Matheval field
Field[2] = MathEval;
Field[2].F = Sprintf("0.0000*F1 +((F1-0e3)/%g)^2 + %g", faultsize, faultsize*10);

//// Equivalent of propagation size on element

// We then define a ‘Threshold’ field, which uses the return value of the
// ‘Distance’ field 1 in order to define a simple change in element size
// depending on the computed distances
//
// SizeMax -                     /------------------
//                              /
//                             /
//                            /
// SizeMin -o----------------/
//          |                |    |
//        Point          DistMin DistMax
//
Field[3] = Threshold;
Field[3].IField = 1;
Field[3].SizeMin = faultsize;
Field[3].SizeMax = boundsize;
Field[3].DistMin = faultsize*10;
Field[3].DistMax = 10e3;
//
Field[4] = Min;
Field[4].FieldsList = {2,3};
Background Field = 2;
