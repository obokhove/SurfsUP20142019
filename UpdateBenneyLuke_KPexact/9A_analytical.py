
import time
import numpy as np
import os.path
from vertical_discr import *
from savings import *
from firedrake import *
from firedrake.petsc import PETSc # 14-04-2020

start_time = time.time()

"""
    ****************************************
    *               Settintgs              *
    **************************************** """
#________________ Kind of data ________________#
#input_data = "coupled"     # coupled with beach
#input_data = "measurements"  # from experiments
input_data = "exact"       # set the wavemaker
#______________ Temporal scheme _______________#
scheme = "SE"
#  "SE": Symplectic-Euler ; "SV": Stormer-Verlet
#__________________ Dimension _________________#
dim = "3D"
#"2D": R(t) and b(x); "3D": R(y,t) and/or b(x,y)
# if input = measurements, the dim must be 2D.
#______ Path and name of the saved files ______#
save_path = 'data/exact_solution/'
fileE = "data/exact_solution/Linf.txt"
outputE = open(fileE,"w")
"""
    *************************************************************
    *                        Saving Files                       *
    ************************************************************* """
save_h = File(os.path.join(save_path, "h.pvd"))

"""
    *************************************************
    *      Definition of the space/time domain      *
    ************************************************* """
if input_data == "exact":
    #______________________ Beach ______________________#
    xb = 7.5                         # Start of the beach
    sb = 0.0                         # Slope of the beach

    #______________________ Basin ______________________#
    g = 9.81                     # Gravitational constant
    H0 = 1.0                # Depth at rest (flat bottom)
    Hend = 1.0              # Depth at the end of the beach
    Lx = 80                   # Length in x 160
    Ly = 40                         # Length in y 80
    Lw = 4.0                     # End of the x-transform

    #____________________ Wavemaker ____________________#
    lamb = 2.0                               # Wavelength
    k = 2*pi/lamb                           # Wave number
    w = sqrt(g*k*tanh(k*H0))             # Wave frequency
    Tw = 2*pi/w                             # Wave period
    gamma = 0.00                         # Wave amplitude
    t_stop = 0.0*Tw          # When to stop the wavemaker

    #______________________ Time _______________________#
    T0 = -3.0                               # Initial time -6.0
    t = T0                            # Temporal variable 6.0
    Tend = 3.0                              # Final time
    dt = 0.1                                # Time-step 0.02

"""
    ****************************************************
    *               Definition of the mesh             *
    **************************************************** """

#_________________ Vertical discretization ________________#
n_z = 8                             # Order of the expansion
Nz = n_z+1                  # Number of points in one element

#________________ Horizontal discretization _______________#
#res = 2.0*10**(-1) #0.05                        # resolution
Nx = 400 + 1               # Number of elements in x
Ny = 400 + 1             # Number of elements in y

#___________________________ Mesh _________________________#
hor_mesh = RectangleMesh(Nx,Ny,Lx,Ly,quadrilateral=True)
hor_mesh.coordinates.dat.data[:,0]+= -0.5*Lx # -0.5*Lx
hor_mesh.coordinates.dat.data[:,1]+= -0.5*Ly-2.0 # -0.5*Ly -2.0 since apparently shifted

"""
    *************************************************
    *       Definition of the function spaces       *
    ************************************************* """
#__________________ For h and psi_1 ___________________#
V = FunctionSpace(hor_mesh, "CG", 1)


"""
    ******************************************************
    *            Definition of the functions             *
    ****************************************************** """
#______________________ At time t^n _____________________#
h_n0 = Function(V)                                   # h^n
                                        # H(x)

"""
    ***********************************************************************************
    *                          Initialisation of the Functions                        *
    ***********************************************************************************"""
#----------------------------------------------------------------------------------------#
#                                       Solutions                                        #
#----------------------------------------------------------------------------------------#
#________________________________________ Depth _________________________________________#

k1 = -1.501
k2 = -0.501
k3 = -0.500
k4 = 0.500
k5 = 0.501
k6 = 1.501

K135 = k1+k3+k5
K235 = k2+k3+k5
K136 = k1+k3+k6
K236 = k2+k3+k6
K145 = k1+k4+k5
K245 = k2+k4+k5
K146 = k1+k4+k6
K246 = k2+k4+k6

KK135 = k1*k1 + k3*k3 + k5*k5
KK235 = k2*k2 + k3*k3 + k5*k5
KK136 = k1*k1 + k3*k3 + k6*k6
KK236 = k2*k2 + k3*k3 + k6*k6
KK145 = k1*k1 + k4*k4 + k5*k5
KK245 = k2*k2 + k4*k4 + k5*k5
KK146 = k1*k1 + k4*k4 + k6*k6
KK246 = k2*k2 + k4*k4 + k6*k6

KKK135 = k1*k1*k1 + k3*k3*k3 + k5*k5*k5
KKK235 = k2*k2*k2 + k3*k3*k3 + k5*k5*k5
KKK136 = k1*k1*k1 + k3*k3*k3 + k6*k6*k6
KKK236 = k2*k2*k2 + k3*k3*k3 + k6*k6*k6
KKK145 = k1*k1*k1 + k4*k4*k4 + k5*k5*k5
KKK245 = k2*k2*k2 + k4*k4*k4 + k5*k5*k5
KKK146 = k1*k1*k1 + k4*k4*k4 + k6*k6*k6
KKK246 = k2*k2*k2 + k4*k4*k4 + k6*k6*k6

T1 = k5*k5*(k3-k1) + k3*k3*(k1-k5) + k1*k1*(k5-k3)
T2 = k5*k5*(k3-k2) + k3*k3*(k2-k5) + k2*k2*(k5-k3)
T3 = k6*k6*(k3-k1) + k3*k3*(k1-k6) + k1*k1*(k6-k3)
T4 = k6*k6*(k3-k2) + k3*k3*(k2-k6) + k2*k2*(k6-k3)
T5 = k5*k5*(k4-k1) + k4*k4*(k1-k5) + k1*k1*(k5-k4)
T6 = k5*k5*(k4-k2) + k4*k4*(k2-k5) + k2*k2*(k5-k4)
T7 = k6*k6*(k4-k1) + k4*k4*(k1-k6) + k1*k1*(k6-k4)
T8 = k6*k6*(k4-k2) + k4*k4*(k2-k6) + k2*k2*(k6-k4)


A_expr = Expression("2.0*( T1*K135*K135*exp( K135*x[0] + KK135*x[1] - KKK135*t) \
                          +T2*K235*K235*exp( K235*x[0] + KK235*x[1] - KKK235*t) \
                          +T3*K136*K136*exp( K136*x[0] + KK136*x[1] - KKK136*t) \
                          +T4*K236*K236*exp( K236*x[0] + KK236*x[1] - KKK236*t) \
                          +T5*K145*K145*exp( K145*x[0] + KK145*x[1] - KKK145*t) \
                          +T6*K245*K245*exp( K245*x[0] + KK245*x[1] - KKK245*t) \
                          +T7*K146*K146*exp( K146*x[0] + KK146*x[1] - KKK146*t) \
                          +T8*K246*K246*exp( K246*x[0] + KK246*x[1] - KKK246*t) )",\
                    K135=K135,K235=K235,K136=K136,K236=K236,K145=K145,K245=K245,K146=K146,K246=K246,\
                    KK135=KK135,KK235=KK235,KK136=KK136,KK236=KK236,KK145=KK145,KK245=KK245,KK146=KK146,KK246=KK246,\
                    KKK135=KKK135,KKK235=KKK235,KKK136=KKK136,KKK236=KKK236,KKK145=KKK145,KKK245=KKK245,KKK146=KKK146,KKK246=KKK246,\
                    T1=T1,T2=T2,T3=T3,T4=T4,T5=T5,T6=T6,T7=T7,T8=T8,t=t)

B_expr = Expression("T1*exp( K135*x[0] + KK135*x[1] - KKK135*t) \
                    +T2*exp( K235*x[0] + KK235*x[1] - KKK235*t) \
                    +T3*exp( K136*x[0] + KK136*x[1] - KKK136*t) \
                    +T4*exp( K236*x[0] + KK236*x[1] - KKK236*t) \
                    +T5*exp( K145*x[0] + KK145*x[1] - KKK145*t) \
                    +T6*exp( K245*x[0] + KK245*x[1] - KKK245*t) \
                    +T7*exp( K146*x[0] + KK146*x[1] - KKK146*t) \
                    +T8*exp( K246*x[0] + KK246*x[1] - KKK246*t) ",\
                    K135=K135,K235=K235,K136=K136,K236=K236,K145=K145,K245=K245,K146=K146,K246=K246,\
                    KK135=KK135,KK235=KK235,KK136=KK136,KK236=KK236,KK145=KK145,KK245=KK245,KK146=KK146,KK246=KK246,\
                    KKK135=KKK135,KKK235=KKK235,KKK136=KKK136,KKK236=KKK236,KKK145=KKK145,KKK245=KKK245,KKK146=KKK146,KKK246=KKK246,\
                    T1=T1,T2=T2,T3=T3,T4=T4,T5=T5,T6=T6,T7=T7,T8=T8,t=t)

dTau_dx_expr = Expression("(T1*K135*exp( K135*x[0] + KK135*x[1] - KKK135*t) \
                          +T2*K235*exp( K235*x[0] + KK235*x[1] - KKK235*t) \
                          +T3*K136*exp( K136*x[0] + KK136*x[1] - KKK136*t) \
                          +T4*K236*exp( K236*x[0] + KK236*x[1] - KKK236*t) \
                          +T5*K145*exp( K145*x[0] + KK145*x[1] - KKK145*t) \
                          +T6*K245*exp( K245*x[0] + KK245*x[1] - KKK245*t) \
                          +T7*K146*exp( K146*x[0] + KK146*x[1] - KKK146*t) \
                          +T8*K246*exp( K246*x[0] + KK246*x[1] - KKK246*t))",\
                          K135=K135,K235=K235,K136=K136,K236=K236,K145=K145,K245=K245,K146=K146,K246=K246,\
                          KK135=KK135,KK235=KK235,KK136=KK136,KK236=KK236,KK145=KK145,KK245=KK245,KK146=KK146,KK246=KK246,\
                          KKK135=KKK135,KKK235=KKK235,KKK136=KKK136,KKK236=KKK236,KKK145=KKK145,KKK245=KKK245,KKK146=KKK146,KKK246=KKK246,\
                          T1=T1,T2=T2,T3=T3,T4=T4,T5=T5,T6=T6,T7=T7,T8=T8,t=t)

dTau_dx = Function(V).interpolate(dTau_dx_expr)

Au = Function(V).interpolate(A_expr)
Bu = Function(V).interpolate(B_expr)
Tau = Function(V).assign(Bu)
Cu = Function(V).assign(2.0*dTau_dx*dTau_dx)

Linf = 0.0
while t<Tend:
    A_expr.t=t
    B_expr.t=t
    dTau_dx_expr.t = t
    dTau_dx.interpolate(dTau_dx_expr)
    Au.interpolate(A_expr)
    Bu.interpolate(B_expr)
    Cu.assign(2.0*dTau_dx*dTau_dx)

    h_n0.interpolate(H0+Au/Bu-Cu/(Bu*Bu))
    # hmmax = max(H0+A_expr/B_expr-2.0*dTau_dx_expr*dTau_dx_expr/(B_expr*B_expr))
    # print(t,hmmax)
    # print(t)
    # 14-04-2020
    with h_n0.dat.vec_ro as hh:
        L_inf_err = hh.max()[1]
        # PETSc.Sys.Print('L_inf error norm = %g' L_inf_err)

    Linf = np.maximum(Linf,L_inf_err)
    print('time = %g, L_inf error norm = %g' % (t, L_inf_err))
    print('%g %g' % (t, L_inf_err),file=outputE)
    
    save_h.write(h_n0)
    t +=dt

print('AbsolutemaxiumLinf %g ' % (Linf),file=outputE)
outputE.close()
print('Linf= %g' % Linf)

