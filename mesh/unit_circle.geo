//+
SetFactory("OpenCASCADE");
Circle(1) = {0, 0, 0, 1, 0, 2*Pi};
//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Physical Surface(1) = {1};  // domain 1

Physical Curve(1) = {1}; // domain 1 boundary

Mesh 2;
RefineMesh;
RefineMesh;
RefineMesh;
RefineMesh;
RefineMesh;

Save "unit_circle.m";

