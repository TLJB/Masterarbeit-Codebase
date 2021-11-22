
// Gmsh project created on Sun Oct 31 20:59:49 2021
startx = -1;
starty = 0;
startz = -1;
h1 = 0.2;
h2 = .5;
width = 2;
height = width;
radius = width /2;
thickness = 0.5;
ih = height/20;
lih = 1;
lh = 4;
lw = 10;
lt = 1;
li = 16;
//+
Point(1) = {startx,starty,startz, h1};
Point(2) = {startx + width/2,starty,startz, h1};
Point(3) = {startx + width,starty,startz, h1};//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Extrude {0, ih, 0} {
  Curve{1}; Curve{2}; Layers {lih}; Recombine;
}
//+
Extrude {0, -ih, 0} {
  Curve{1}; Curve{2}; Layers {lih}; Recombine;
}
//+
Extrude {0, height/2 - ih, 0} {
  Curve{3}; Curve{7}; Layers {lh}; Recombine;
}
//+
Extrude {0, -height/2 + ih, 0} {
  Curve{11}; Curve{15}; Layers {lh}; Recombine;
}
//+
Physical Curve("interface", 35) = {2};
//+
Physical Surface("bulk_2", 36) = {6, 10, 26, 22};
//+
Physical Surface("bulk_1", 37) = {30, 34, 18, 14};
//+
Transfinite Curve {1, 2} = li Using Progression 1;
