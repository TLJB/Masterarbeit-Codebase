// Gmsh project created on Sun Oct 31 20:59:49 2021
startx = -1;
starty = -1;
startz = -1;
h1 = .5;
h2 = 2;
height = 2;
width = 2;
radius = width /2;
thickness = 2;
lh = 10;
lw = 10;
lt = 10;
li = 12;
//+
Point(1) = {startx,starty,startz, h1};
Point(2) = {startx + width, starty,startz, h2};
Point(3) = {startx, starty+radius,startz, h1};
Point(4) = {startx + radius, starty,startz, h1};
Point(5) = {startx + width, starty+height,startz, h2};
Point(6) = {startx, starty+height,startz, h2};//+
Point(7) = {startx, starty,startz + radius, h1};

Point(8) = {startx,starty,startz + thickness, h1};
Point(9) = {startx + width, starty,startz + thickness, h2};
Point(10) = {startx + width, starty+height,startz + thickness, h2};
Point(11) = {startx, starty+height,startz + thickness, h2};//+
Line(1) = {1, 4};
//+
Line(2) = {4, 2};
//+
Line(3) = {2, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {6, 3};
//+
Line(6) = {3, 1};
//+
Circle(7) = {3, 1, 4};
//+
Circle(8) = {4, 1, 7};
//+
Circle(9) = {3, 1, 7};
//+
Line(10) = {1, 7};
//+
Curve Loop(1) = {7, -1, -6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, -10, -6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, 8, -10};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {7, 8, -9};
//+
Surface(4) = {4};
//+
Surface Loop(1) = {3, 1, 4, 2};
//+
Volume(1) = {1};
//+
Physical Volume("bulk_2", 11) = {1};
//+
Line(11) = {8, 9};
//+
Line(12) = {9, 10};
//+
Line(13) = {10, 11};
//+
Line(14) = {11, 8};
//+
Line(15) = {7, 8};
//+
Line(16) = {6, 11};
//+
Line(17) = {2, 9};
//+
Line(18) = {5, 10};
//+
Curve Loop(5) = {2, 17, -11, -15, -8};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {16, -13, -18, 4};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {5, 9, 15, -14, -16};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {14, 11, 12, 13};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {12, -18, -3, 17};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {4, 5, 7, 2, 3};
//+
Plane Surface(10) = {10};
//+
Surface Loop(2) = {6, 7, 10, 5, 9, 8, 4};
//+
Volume(2) = {2};
//+
Physical Surface("interface", 19) = {4};
//+
Physical Volume("bulk_1", 20) = {2};
//+
Transfinite Curve {7, 8, 9} = li Using Progression 1;
