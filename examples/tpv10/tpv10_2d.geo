SetFactory("OpenCASCADE");
faultsize = 200;
boundsize = 5000;

R = 75e3; // semi-circle as remote boundary
FaultLen = 15e3;
FaultDipAngle = 60;
FaultDipAngle = FaultDipAngle/180.0*Pi;
Point(1) = {0, 0, 0, 1.0};
Point(2) = {R, 0, 0, 1.0};
Point(3) = {-R, 0, 0, 1.0};
Circle(5) = {2, 1, 3};
Point(4) = {FaultLen*Cos(FaultDipAngle), -FaultLen*Sin(FaultDipAngle), 0, 1.0};
Line(1) = {1, 2};
Line(2) = {1, 4};
Line(3) = {1, 3};

MeshSize {:} = boundsize;
MeshSize {1, 4} = faultsize;
Curve Loop(1) = {3, -5, -1};
Surface(1) = {1};
//
BooleanFragments{ Surface{1}; Delete; }{ Curve{2}; Delete; }

//Mesh.Algorithm = 6; // Frontal-Delaunay
Mesh.Algorithm = 9; //  Packing of Parallelograms
Mesh.RecombinationAlgorithm = 2;
Mesh.RecombineAll = 1;
//+

