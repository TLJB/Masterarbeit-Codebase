// Gmsh project created on Sun Oct 31 20:59:49 2021
startx = -1;
starty = -1;
startz = -1;
h1 = .5;
h2 = 1;
height = 2;
width = 2;
radius = width /4;
thickness = 2;
middlex = startx + width/2;
middley = starty + height/2;
middlez = startz + 0*thickness/2;

lh = 10;
lw = 4;
lt = 10;
li = 12;
//+
Point(1) = {startx,starty,startz, h1};
Point(2) = {startx + width, starty,startz, h2};
Point(3) = {startx + width, starty+height,startz, h2};
Point(4) = {startx, starty+height,startz, h2};//+

Point(5) = {middlex, middley, middlez, h2};//+
Point(6) = {middlex - radius, middley, middlez, h2};//+
Point(7) = {middlex + radius, middley, middlez, h2};//+
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {6, 5, 7};
//+
Circle(6) = {7, 5, 6};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Curve Loop(2) = {6, 5};
//+
Plane Surface(1) = {1, 2};
//+
Plane Surface(2) = {2};
//+
Physical Curve("interface", 7) = {6, 5};
//+
Physical Surface("bulk_1", 8) = {1};
//+
Physical Surface("bulk_2", 9) = {2};
//+
Transfinite Curve {6, 5} = li Using Progression 1;
//+
Transfinite Curve {4, 3, 2, 1} = lw Using Progression 1;
