
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
Point(2) = {startx + width,starty,startz, h1};//+
Line(1) = {1, 2};
//+
Extrude {0, ih, 0} {
  Curve{1}; Layers {lih}; Recombine;
}
//+
Extrude {0, ih, 0} {
  Curve{1}; Layers {lih}; Recombine;
}
//+
Extrude {0, -ih, 0} {
  Curve{1}; Layers {lih}; Recombine;
}
//+
Extrude {0, height/2 - ih, 0} {
  Curve{2}; Layers {lh}; Recombine;
}
//+
Extrude {0, -height/2 + ih, 0} {
  Curve{6}; Layers {lh}; Recombine;
}//+
//+
Physical Curve("interface", 18) = {1};
//+
Physical Surface("bulk_1", 19) = {17, 9};
//+
Physical Surface("bulk_2", 20) = {5, 13};
//+
Transfinite Curve {1} = li Using Progression 1;
