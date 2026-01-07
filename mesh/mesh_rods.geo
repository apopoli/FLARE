// Parametri
L = 0.1; // larghezza
H = 0.02; // altezza
d = 0.006; // distanza tra gli elettrodi

r = 0.0015; // raggio di arrotondamento (0 < r < min(L/2, H))
lms = 0.001/2; // lunghezza mesh

// centrare tutto nell'origine
x0 = -(L+d/2); // origine x
y0 = 0; // origine y

p1 = newp;
Point(p1) = {x0, y0+H, 0, lms*3};
p2 = newp;
Point(p2) = {x0, y0, 0, lms*3};

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

// centro circonf SE
xp1 = x0+L; yp1 = y0+0+r;
pc21 = newp;
Point(pc21) = {xp1, yp1, 0, lms};
xp2 = x0+L-r;
yp2 = y0+0;
pc22 = newp;
Point(pc22) = {xp2, yp2, 0, lms};
// centro circonf SE
pc20 = newp;
Point(pc20) = {xp1-r,yp2+r, 0, lms};
Circle(2) = {pc21, pc20, pc22};

//+
Line(3) = {p1, pc11};
//+
Line(4) = {pc12, pc21};
//+
Line(5) = {pc22, p2};
//+
Line(6) = {p2, p1};
//+
Curve Loop(1) = {3, -1, 4, 2, 5, 6};
//+
Plane Surface(1) = {1};
//+

Symmetry {1, 0, 0, -d/4} {
  Duplicata { Point{3}; Point{4}; Point{5}; Point{6}; Point{8}; Point{7}; Point{1}; Point{2}; Curve{1}; Curve{4}; Curve{2}; Curve{3}; Curve{5}; Curve{6}; Surface{1}; }
}
