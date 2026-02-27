// Parametri
L = 0.07; // larghezza
H = 0.01; // altezza
H_edge = 3*H; // altezza bordo superiore
d = 0.006; // distanza tra gli elettrodi

If (!Exists(r))
  r = 1.5e-4; // raggio di arrotondamento (0 < r < min(L/2, H))
EndIf

lms = 0.001/1.5; // lunghezza mesh

// centrare tutto nell'origine
x0 = -(L+d/2); // origine x
y0 = 0; // origine y

// centrare tutto nell'origine
x0 = -(L+d/2); // origine x
y0 = 0; // origine y

p0 = newp;
Point(p0) = {x0, y0+H_edge, 0, lms*5};

p1 = newp;
Point(p1) = {x0, y0+H, 0, lms*3};

xp1 = x0+L-r; yp1 = y0+H;
pc11 = newp;
Point(pc11) = {xp1, yp1, 0, lms};
xp2 = x0+L; yp2 = y0+H-r;
pc12 = newp;
Point(pc12) = {xp2, yp2, 0, lms};

// centro circonf NE
pc10 = newp;
Point(pc10) = {xp2-r,yp1-r, 0, lms};
Circle(1) = {pc12, pc10, pc11};

p2 = newp;
Point(p2) = {x0+L, y0, 0, lms};

Symmetry {1, 0, 0, 0} {
  Duplicata { Point{p0}; Point{p1}; Point{pc11}; Point{pc12}; Point{pc10}; Point{p2}; Curve{1}; }
}
//+
Line(3) = {1, 2};
//+
Line(4) = {2, 3};
//+
Line(5) = {4, 6};
//+
Line(6) = {6, 12};
//+
Line(7) = {12, 10};
//+
Line(8) = {9, 8};
//+
Line(9) = {8, 7};
//+
Line(10) = {7, 1};
//+
Curve Loop(1) = {10, 3, 4, -1, 5, 6, 7, 2, 8, 9};
//+
Plane Surface(1) = {1};
//+
Physical Curve(11) = {4, 1, 5};
//+
Physical Curve(12) = {7, 2, 8};
//+
Physical Curve(13) = {3, 10, 9, 6};

Physical Surface(1) = {1};

//+
Transfinite Curve {1,2} = 25 Using Progression 1;
//+
//+
Transfinite Curve {5, -7} = 60 Using Progression 1.05;

Transfinite Curve {6} = 20 Using Progression 1;

Transfinite Curve {4} = 95 Using Progression 0.95;
Transfinite Curve {8} = 95 Using Progression 1.05;

// MESH 2D
Mesh 2;
Mesh.MshFileVersion = 2;

Save "mesh_rods_half.m";

