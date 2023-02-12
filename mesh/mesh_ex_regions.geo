// Gmsh project created on Fri Sep 28 18:07:57 2018
SetFactory("OpenCASCADE");
//dimension of external box
a = 2;
// circle 1 center and radius
xc1 = 0.1;
yc1 = 0.2;
R1 = 0.04;
// circle 2 center and radius
xc2 = -0.5;
yc2 = -0.4;
R2 = 0.04;
//box SW corner coordinates, Lx & Ly 
xb = -0.1;
yb = -0.2;
Lx = 0.5;
Ly = 0.05;
//+ box vertices definition
Point(1) = {-a, -a, 0, 0.2};
//+
Point(2) = {a, -a, 0, 0.2};
//+
Point(3) = {a, a, 0, 0.2};
//+
Point(4) = {-a, a, 0, 0.2};
//
//+ circle 1 points and center definition
Point(5) = {xc1+R1, yc1, 0, .05};
//+
Point(6) = {xc1, yc1+R1, 0, .05};
//+
Point(7) = {xc1-R1, yc1, 0, .05};
//+
Point(8) = {xc1, yc1-R1, 0, .05};
//+
Point(9) = {xc1, yc1, 0, .1};
//
//+ circle 2 points and center definition
Point(10) = {xc2+R2, yc2, 0, .05};
//+
Point(11) = {xc2, yc2+R2, 0, .05};
//+
Point(12) = {xc2-R2, yc2, 0, .05};
//+
Point(13) = {xc2, yc2-R2, 0, .05};
//+
Point(14) = {xc2, yc2, 0, .1};
//+ box definition
Point(15) = {xb, yb, 0, .005};
//+
Point(16) = {xb+Lx, yb, 0, .005};
//+
Point(17) = {xb+Lx, yb+Ly, 0, .005};
//+
Point(18) = {xb, yb+Ly, 0, .005};
//
//+ external box edges
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//
//+ circle 1 arcs
Circle(5) = {5, 9, 6};
//+
Circle(6) = {6, 9, 7};
//+
Circle(7) = {7, 9, 8};
//+
Circle(8) = {8, 9, 5};
//
//+ circle 2 arcs
Circle(9) = {10, 14, 11};
//+
Circle(10) = {11, 14, 12};
//+
Circle(11) = {12, 14, 13};
//+
Circle(12) = {13, 14, 10};
//
//+ box edges
Line(13) = {15, 16};
//+
Line(14) = {16, 17};
//+
Line(15) = {17, 18};
//+
Line(16) = {18, 15};
//
//+ external boundary closed curve definition
Curve Loop(1) = {4, 1, 2, 3};
//
//+ cirle 1 closed curve definition
Curve Loop(2) = {5, 6, 7, 8};
//
//+ cirle 2 closed curve definition
Curve Loop(3) = {9, 10, 11, 12};
//
//+ box closed curve definition
Curve Loop(4) = {13, 14, 15, 16};
//+ 
Plane Surface(1) = {1, 2, 3, 4};
//+
Plane Surface(2) = {2};
//+
Plane Surface(3) = {3};
//+
Plane Surface(4) = {4};
//+
Transfinite Curve {6, 7, 8, 5} = 15 Using Progression 1;
//+
Transfinite Curve {9, 10, 11, 12} = 15 Using Progression 1;
//+
Transfinite Curve {13} = 61 Using Bump .2;
//+
Transfinite Curve {14} = 32 Using Bump .333;
//+
Transfinite Curve {15} = 61 Using Bump .2;
//+
Transfinite Curve {16} = 32 Using Bump .333;
//+
Transfinite Curve {1, 2, 3, 4} = 21 Using Progression 1;
//+
Transfinite Surface {4} = {15, 16, 17, 18};
//+
Physical Surface("medium") = {1};
//+
Physical Surface("circle1") = {2};
//+
Physical Surface("circle2") = {3};
//+
Physical Surface("box") = {4};
//+
Physical Curve("S") = {1};
//+
Physical Curve("E") = {2};
//+
Physical Curve("N") = {3};
//+
Physical Curve("W") = {4};

Mesh 2;
Mesh.MshFileVersion = 2;

Save "mesh_ex_regions.m";
