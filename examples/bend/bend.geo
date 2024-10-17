// Gmsh project created on Wed Sep 22 12:42:09 2021

SetFactory("OpenCASCADE");
faultsize = 0.2e3;
boundsize = 10e3;
R = 100e3;
Point(1) = {20e3, 0, 0, 1.0};
BendAngle = -10;
//BendAngle = 5;
//BendAngle = 10;
BendAngle = BendAngle/180*Pi;
FaultLen2 = 20e3;
Point(2) = {20e3+FaultLen2*Cos(BendAngle), FaultLen2*Sin(BendAngle), 0, 1.0};
Point(3) = {-40e3, 0, 0, 1.0};
Line(1) = {3, 1};
Line(2) = {1, 2};
Circle(4) = {0, 0, 0, R, 0, 2*Pi};
Curve Loop(1) = {4};
Surface(1) = {1};
BooleanFragments{ Surface{1}; Delete; }{ Curve{1:2}; Delete; }
//+
MeshSize {:} = boundsize;
MeshSize {1:3} = faultsize;
//+
Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.Algorithm = 9; //  Packing of Parallelograms
Mesh.RecombinationAlgorithm = 2;
Mesh.RecombineAll = 1;

Field[1] = Distance;
Field[1].CurvesList = {1:2};

// Matheval field
Field[2] = MathEval;
Field[2].F = Sprintf("0.0000*F1 +((F1-0e3)/%g)^2 + %g", faultsize, faultsize);

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
Field[3].DistMax = R-boundsize*4;
//
Field[4] = Min;
Field[4].FieldsList = {2,3};
Background Field = 3;
