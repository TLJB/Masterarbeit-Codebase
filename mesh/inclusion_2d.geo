// Gmsh project created on Sun Oct 31 20:59:49 2021
startx = 0;
starty = 0;
startz = 0;
h1 = .5;
h2 = 1;
height = 2;
width = 2;
radius = width /2;
thickness = 2;
lh = 10;
lw = 2;
lt = 10;
li = 12;
//+
Point(1) = {startx,starty,startz, h1};
Point(2) = {startx + width, starty,startz, h2};
Point(3) = {startx, starty+radius,startz, h1};
Point(4) = {startx + radius, starty,startz, h1};
Point(5) = {startx + width, starty+height,startz, h2};
Point(6) = {startx, starty+height,startz, h2};//+

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
Circle(7) = {3, 1, 4};//+
Curve Loop(1) = {7, -1, -6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, 2, 3, 4, 5};
//+
Plane Surface(2) = {2};
//+
Physical Curve("interface", 8) = {7};
//+
Physical Surface("bulk_1", 9) = {1};
//+
Physical Surface("bulk_2", 10) = {2};
//+
Transfinite Curve {7} = li Using Progression 1;
//+
// delete for quarter//+
// Symmetry {1, 0, 0, 0} {
//   Duplicata { Surface{2}; Surface{1}; }
// }
// //+
// Symmetry {0, 1, 0, 0} {
//   Duplicata { Surface{8}; Surface{14}; Surface{1}; Surface{2}; }
// }
// //+
// //+
// Transfinite Curve {7, 28, 18, 9} = li Using Progression 1;
// //+
// Transfinite Curve {12, 4, 3, 34, 35, 21, 20, 11} = lw Using Progression 1;
// //+
// Physical Curve("interface", 36) = {7, 28, 18, 9};
// //+
// Physical Surface("bulk_2", 37) = {2, 17, 14, 27};
// //+
// Physical Surface("bulk_1", 38) = {8, 31, 23, 1};
