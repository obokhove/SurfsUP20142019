from firedrake import *

n = 10
nlayers = 10
H0 = 1.0
m = UnitIntervalMesh(n)
mesh = ExtrudedMesh(m, nlayers, layer_height=H0/nlayers)

V = FunctionSpace(mesh, "CG", 1)
S = FunctionSpace(mesh, "CG", 1, vfamily="R", vdegree=0)

# zbar
Wzb = MixedFunctionSpace((V))
wzb = Function(Wzb)
a = Constant(2.0)

z, zbar, lamda, lamda0 = split(w)

S = ( zbar.dx(1)*zbar.dx(1)/2. )*dx # test simple 2nd order case first

# S = (z*z + z.dx(1)*z.dx(1) +
#     zbar*zbar + zbar.dx(1)*zbar.dx(1) +
#     lamda*lamda + lamda0*lamda0)*dx + a*lamda0*ds_t

zbar = TrialFunctions(Wzb)

#Dirichlet condition
bcs = [DirichletBC(Wzb.sub(0), 0., "bottom")]

eqn = derivative(S, wzb)

prob = NonlinearVariationalProblem(eqn, wzb, bcs=bcs)

fieldsplit_parameters= {"ksp_type":"gmres",
                        "pc_type": "fieldsplit",
                        "pc_fieldsplit_type": "additive",
                        "ksp_monitor":None,
                        "ksp_converged_reason":None,
                        "fieldsplit_0_ksp_type":"preonly",
                        "fieldsplit_0_pc_type":"lu"}

sor_parameters = {"ksp_type":"gmres",
                  "pc_type": "sor",
                  "ksp_monitor":None,
                  "mat_type":"aij",
                  "ksp_converged_reason":None}   

solver = NonlinearVariationalSolver(prob,solver_parameters=fieldsplit_parameters)

solver.solve()

zbar = wzb.split()

File('isolation.pvd').write(zbar)
