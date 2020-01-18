x0 = 0;
z0 = 0;

L = 1;
H = 1;
H0 = 0.5;

Nx = 40;
Nz = 80;

lc = 1e-2;

prog = 1.07; // progression

p1 = newp; Point(p1) = {x0, z0, 0, lc};
p2 = newp; Point(p2) = {x0+L, z0, 0, lc};
p3 = newp; Point(p3) = {x0+L, z0+H0, 0, lc};
p4 = newp; Point(p4) = {x0, z0+H0, 0, lc};

p5 = newp; Point(p5) = {x0+L, z0+H, 0, lc};
p6 = newp; Point(p6) = {x0, z0+H, 0, lc};

l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p3,p2};
l3 = newl; Line(l3) = {p3,p4};
l4 = newl; Line(l4) = {p4,p1};

l5 = newl; Line(l5) = {p3,p5};
l6 = newl; Line(l6) = {p5,p6};
l7 = newl; Line(l7) = {p4,p6};

// Create surface
loop1 = newreg; Line Loop(loop1) = {l1,-l2,l3,l4};
s1 = newreg;
Surface(s1) = {loop1};

loop2 = newreg; Line Loop(loop2) = {l3,l7,-l6,-l5};
s2 = newreg;
Surface(s2) = {loop2};

Transfinite Line{l1} = Nx+1;
Transfinite Line{l3} = Nx+1;
Transfinite Line{l2} = Nz+1 Using Progression prog;
Transfinite Line{l4} = Nz+1 Using Progression prog;

Transfinite Line{l6} = Nx+1;
Transfinite Line{l5} = Nz+1 Using Progression prog;
Transfinite Line{l7} = Nz+1 Using Progression prog;

Transfinite Surface{s1} = {p1,p2,p3,p4};
Transfinite Surface{s2} = {p3,p5,p6,p4};
Recombine Surface{s1};
Recombine Surface{s2};

// Apply an elliptic smoother to the grid
//Mesh.Smoothing = 100;

Physical Surface(1) = s1;
Physical Surface(2) = s2;
