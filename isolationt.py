from firedrake import *

n = 10
nlayers = 10
H0 = 1.0
m = UnitIntervalMesh(n)
mesh = ExtrudedMesh(m, nlayers, layer_height=H0/nlayers)

V = FunctionSpace(mesh, "CG", 1)
S = FunctionSpace(mesh, "CG", 1, vfamily="R", vdegree=0)

# z, zbar, lambda, lambda0
W = MixedFunctionSpace((V,V,V,S))
w = Function(W)
a = Constant(2.0)

z, zbar, lamda, lamda0 = split(w)

S = (zbar.dx(1)*zbar.dx(1)/2. +
     lamda*(z-zbar))*dx + lamda0*(z-a)*ds_t

#S = (z*z + z.dx(1)*z.dx(1) +
#     zbar*zbar + zbar.dx(1)*zbar.dx(1) +
#     lamda*lamda + lamda0*lamda0)*dx + a*lamda0*ds_t

z, zbar, lamda, lamda0 = TrialFunctions(W)
dz, dzbar, dlamda, dlamda0 = TestFunctions(W)

Jp = (
    z*dz + z.dx(1)*dz.dx(1) +
    zbar*dzbar + zbar.dx(1)*dzbar.dx(1) +
    lamda*dlamda
)*dx + lamda0*dlamda0*dx

#Dirichlet condition
bcs = [DirichletBC(W.sub(1), 0., "bottom")]

eqn = derivative(S, w)

prob = NonlinearVariationalProblem(eqn, w, bcs=bcs, Jp=Jp)

fieldsplit_parameters= {"ksp_type":"gmres",
                        "pc_type": "fieldsplit",
                        "pc_fieldsplit_type": "additive",
                        "ksp_monitor":None,
                        "ksp_converged_reason":None,
                        "fieldsplit_0_ksp_type":"preonly",
                        "fieldsplit_0_pc_type":"lu",
                        "fieldsplit_1_ksp_type":"preonly",
                        "fieldsplit_1_pc_type":"lu",
                        "fieldsplit_2_ksp_type":"preonly",
                        "fieldsplit_2_pc_type":"lu",
                        "fieldsplit_3_ksp_type":"preonly",
                        "fieldsplit_3_pc_type":"lu"}

sor_parameters = {"ksp_type":"gmres",
                  "pc_type": "sor",
                  "ksp_monitor":None,
                  "mat_type":"aij",
                  "ksp_converged_reason":None}   

solver = NonlinearVariationalSolver(prob,solver_parameters=
                                    fieldsplit_parameters)

solver.solve()

z, zbar, lamda, lamda0 = w.split()

File('isolation.pvd').write(z, zbar, lamda)
