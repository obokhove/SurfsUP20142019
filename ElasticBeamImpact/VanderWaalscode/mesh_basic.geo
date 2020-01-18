x0 = 0;
y0 = 0;
z0 = 0;

L = 1;
H = 1;

Nx = 10;
Nz = 20;

lc = 1e-2;

p1 = newp; Point(p1) = {x0, y0, z0, lc};
p2 = newp; Point(p2) = {x0+L, y0, z0, lc};
p3 = newp; Point(p3) = {x0+L, y0, z0+H, lc};
p4 = newp; Point(p4) = {x0, y0, z0+H, lc};

l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p3,p2};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p4,p1};

// Create surface
loop1 = newreg; Line Loop(loop1) = {l1,-l2,l3,l4};
s1 = newreg;
Surface(s1) = {loop1};


Transfinite Line{l1} = Nx+1;
Transfinite Line{l3} = Nx+1;
Transfinite Line{l2} = Nz+1;
Transfinite Line{l4} = Nz+1;


Transfinite Surface{s1} = {p1,p2,p3,p4};
//Recombine Surface{s1};

// Apply an elliptic smoother to the grid
//Mesh.Smoothing = 100;

Physical Surface(1) = s1;
