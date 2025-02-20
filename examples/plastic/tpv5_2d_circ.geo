//+
SetFactory("OpenCASCADE");
faultsize = 500;
boundsize = 20e3;
R = 100e3;
Point(1) = {-40e3, 0, 0, 1.0};
Point(2) = {+40e3, 0, 0, 1.0};
Line(1) = {1, 2};
Circle(4) = {0, 0, 0, R, 0, 2*Pi};
Curve Loop(1) = {4};
Surface(1) = {1};
//+
BooleanFragments{ Surface{1}; Delete; }{ Curve{1}; Delete; }

MeshSize {:} = boundsize;
MeshSize {1,2} = faultsize;

// "Frontal-Delaunay" 2D meshing algorithm (Mesh.Algorithm = 6)
// usually leads to the highest quality meshes, the
// "Delaunay" algorithm (Mesh.Algorithm = 5) will handle complex mesh size


Mesh.Algorithm = 5; // Frontal-Delaunay
//Mesh.Algorithm = 9; //  Packing of Parallelograms
// Mesh recombination algorithm (0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad)
Mesh.RecombinationAlgorithm = 2;
//Mesh.RecombineAll = 1;

Mesh.SubdivisionAlgorithm = 1;
//Mesh subdivision algorithm (0: none, 1: all quadrangles, 2: all hexahedra, 3:barycentric)
//RefineMesh;
//Physical Surface("fault") = {1};

//
Field[1] = Distance;
Field[1].CurvesList = {1};

// Matheval field
Field[2] = MathEval;
Field[2].F = Sprintf("0.0000*F1 +((F1-0e3)/0.3e3)^2 + %g", faultsize);

//// Equivalent of propagation size on element
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
Field[3].DistMin = faultsize*5;
Field[3].DistMax = R-boundsize*1;
Field[3].DistMax = R-10e3;
//
Field[4] = Min;
Field[4].FieldsList = {2,3};
Background Field = 3;
