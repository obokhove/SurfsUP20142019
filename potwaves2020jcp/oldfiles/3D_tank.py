
import time
import numpy as np
import os.path
from vertical_discr import *
from savings import *
from Settings import *
from firedrake import *
import solvers as DW_solvers

start_time = time.time()


"""
    ****************************************
    *               Settings               *
    **************************************** """
input_data, scheme, dim, save_path = test_case()
H0, xb, sb, H_expr, Hend, Lx, Ly, Lw, res_x, res_y, n_z = domain()
T0, t, dt, Tend, dt_save = set_time()
g, lamb, k, w, Tw, gamma, t_stop, WM_expr, dWM_dt_expr, dWM_dy_expr = wavemaker(dim, H0, Ly, Lw, t)

print('Settings loaded')
print('Creation of the mesh...')

"""
    ****************************************************
    *               Definition of the mesh             *
    **************************************************** """

#_________________ Vertical discretization ________________#
Nz = n_z+1                  # Number of point in one element

#________________ Horizontal discretization _______________#
Nx = round(Lx/res_x)                 # Number of elements in x
Ny = round(Ly/res_y)                 # Number of elements in y

#___________________________ Mesh _________________________#
if dim=="2D":                                   #(x,z)-waves
    hor_mesh = IntervalMesh(Nx,Lx)
else:                                         #(x,y,z)-waves
    hor_mesh = RectangleMesh(Nx,Ny,Lx,Ly,quadrilateral=True)

print('Mesh created...')
print('Definition of the functions...')
"""
    *************************************************
    *       Definition of the function spaces       *
    ************************************************* """
#___________________ For h and psi_1 ___________________#
V = FunctionSpace(hor_mesh, "CG", 1)
#_____________________ For hat_psi _____________________#
Vec = VectorFunctionSpace(hor_mesh, "CG", 1,dim=n_z)
#_________________ Mixed function space ________________#
V_mixed = V*Vec # to solve simultaneous weak formulations



"""
    ******************************************************
    *            Definition of the functions             *
    ****************************************************** """

if scheme=="SE": #_________ Symplectic-Euler scheme _________#
    #______________________ At time t^n _____________________#
    h_n0 = Function(V)                                   # h^n
    psi_1_n0 = Function(V)                           # psi_1^n
    hat_psi_n0 = Function(Vec)                     # hat_psi^n

    #________________ At time t^{n+1} and t^* _______________#
    psi_1_n1 = Function(V)                       # psi_1^{n+1}
    w_n1 = Function(V_mixed)
    h_n1, hat_psi_star = split(w_n1)      # h^{n+1}, hat_psi^*
    hat_psi_n1 = Function(Vec)    # to visualise hat_psi^{n+1}
else: #________________ Stormer-Verlet scheme _______________#
    #______________________ At time t^n _____________________#
    h_n0 = Function(V)                                   # h^n
    psi_1_n0 = Function(V)                           # psi_1^n
    hat_psi_n0 = Function(Vec)                     # hat_psi^n

    #_______________ At time t^{n+1/2} and t^* ______________#
    w_half = Function(V_mixed)        # to obtain psi^{n+1/2},
    psi_1_half, hat_psi_star = split(w_half)   # and hat_psi^*

    #_______________ At time t^{n+1} and t^** _______________#
    psi_1_n1 = Function(V)                       # psi_1^{n+1}
    w_n1 = Function(V_mixed)              # to obtain h^{n+1},
    h_n1, hat_psi_aux = split(w_n1)         # and hat_psi^{**}
    hat_psi_n1 = Function(Vec)    # to visualise hat_psi^{n+1}


#_______________________ x coordinate _______________________#
x_coord = Function(V).interpolate(Expression('x[0]'))

#___________________________ Beach __________________________#
beach = Function(V)                                     # b(x)

#_______________________ Depth at rest ______________________#
H = Function(V)                                         # H(x)

#_________________________ Wavemaker ________________________#
WM = Function(V)                                  # R(x,y;t^n)
dWM_dt = Function(V)                               # (dR/dt)^n
dWM_dy = Function(V)                               # (dR/dy)^n
if scheme=="SV":                         # For Stormer-Verlet:
    WM_half = Function(V)                   # R(x,y;t^{n+1/2})
    dWM_half_dt = Function(V)                # (dR/dt)^{n+1/2}
    dWM_half_dy = Function(V)                # (dR/dy)^{n+1/2}
WM_n1 = Function(V)                           # R(x,y;t^{n+1})
dWM_n1_dt = Function(V)                        # (dR/dt)^{n+1}
dWM_n1_dy = Function(V)                        # (dR/dy)^{n+1}

#______________________ Trial functions _____________________#
psi_1 = TrialFunction(V)      # psi_1^{n+1} for linear solvers
hat_psi = TrialFunction(Vec)# hat_psi^{n+1} for linear solvers

#_______________________ Test functions _____________________#
delta_h = TestFunction(V)                         # from dH/dh
delta_hat_psi = TestFunction(Vec)           # from dH/dhat_psi
w_t = TestFunction(V_mixed)                # from dH/dpsi_1...
delta_psi, delta_hat_star = split(w_t)    # ...and dH/dhat_psi

print('...functions created')
print(' Initalisation of the functions...')
"""
    ***********************************************************************************
    *                          Initialisation of the Functions                        *
    ***********************************************************************************"""
#---------------------------- Topography ----------------------------#
H.interpolate(H_expr)                             # Depth at rest H(x)
beach.interpolate(H0-H)                                   # Beach b(x)
#----------------------------------------------------------------------------------------#
#                                       Wavemaker                                        #
#----------------------------------------------------------------------------------------#
if input_data=="measurements":#------------- Interpolate measurements -------------#
    #------------------------------- Wavemaker motion -------------------------------#
    WM_expr = \
    Expression("((wm2*(t-t1) - wm1*(t-t2))/(t2-t1))*0.5*(1.0+copysign(1.0,Lw-x[0]))",\
               wm2=wm_data[1], wm1=wm_data[0], t1=t_data[0], t2=t_data[1], t=t, Lw=Lw)
    WM.interpolate(WM_expr)

    #------------------------------- Wavemaker velocity -----------------------------#
    dWM_dt_expr = \
    Expression("((dwm2*(t-t1) - dwm1*(t-t2))/(t2-t1))*0.5*(1.0+copysign(1.0,Lw-x[0]))",\
               dwm2=wm_vel_data[1], dwm1=wm_vel_data[0], t1=t_data[0], t2=t_data[1], \
               t=t, Lw=Lw)
    dWM_dt.interpolate(dWM_dt_expr)
else:
    WM.interpolate(WM_expr)                  # \tilde{R}(x,y;t)
    dWM_dt.interpolate(dWM_dt_expr)          # d\tilde{R}/dt
    dWM_dy.interpolate(dWM_dy_expr)          # d\tilde{R}/dy


#----------------------------------------------------------------------------------------#
#                                       Solutions                                        #
#----------------------------------------------------------------------------------------#
#________________________________________ Depth _________________________________________#
h_n0.assign(H)                                                         # h(x,y;t=0) = H(x)
w_n1.sub(0).assign(H) # h^{n+1}

#_____________________ Velocity pot. at the surface: phi(x,y,z=h;t) _____________________#
psi_1_n0.assign(0.0)                                                 # \psi_1(x,y;t=0) = 0

#_____________________ Velocity pot. in depth: phi(x,y,z<h;t) _____________________#
for i in range(0,n_z):
    hat_psi_n0.dat.data[:,i] = 0.0 # psi_i^n
    w_n1.sub(1).dat.data[:,i] = 0.0  # psi_i^{*}


if scheme=="SV": #______________________________ Stormer-Verlet scheme _____________________________#
    w_half.sub(0).assign(0.0)
    for i in range(0,n_z):
        w_half.sub(1).dat.data[:,i]= 0.0 # whalf is psi_1 and hat_psi^* for SV (not used for SE)



print('...functions initialised')
print('Assembling z-matrices...')
"""
    ************************
    * Compute the matrices *
    ************************ """
#_______ Initialization ______#
A = np.eye(Nz,Nz)*0.0
M = np.eye(Nz,Nz)*0.0
D = np.eye(Nz,Nz)*0.0
S = np.eye(Nz,Nz)*0.0
Ik = np.eye(Nz,1)*0.0

#____ Filling the matrices ___#
for i in range(0,Nz):
    for j in range(0,Nz):
        A[i,j]=A_ij(i,j,n_z,H0)
        M[i,j]=M_ij(i,j,n_z,H0)
        D[i,j]=D_ij(i,j,n_z,H0)
        S[i,j]=S_ij(i,j,n_z,H0)
    Ik[i] = I_i(i,n_z,H0)

#________ Submatrices ________#
A11 = A[0,0]
A1N = as_tensor(A[0,1:])
AN1 = as_tensor(A[1:,0])
ANN = as_tensor(A[1:,1:])

M11 = M[0,0]
M1N = as_tensor(M[0,1:])
MN1 = as_tensor(M[1:,0])
MNN = as_tensor(M[1:,1:])

D11 = D[0,0]
D1N = as_tensor(D[0,1:])
DN1 = as_tensor(D[1:,0])
DNN = as_tensor(D[1:,1:])

S11 = S[0,0]
S1N = as_tensor(S[0,1:])
SN1 = as_tensor(S[1:,0])
SNN = as_tensor(S[1:,1:])

I1 = Ik[0,0]
IN=as_tensor(Ik[1:,0])

print('... end of assembling')
print('Initialisation of the solvers...')
"""
    ************************************************************************************************************************
    *                                                   Weak Formulations                                                  *
    ************************************************************************************************************************ """

if scheme=="SE": #_____________________________________________ Symplectic-Euler ______________________________________________#
    #------------------------ Step 1 : Update h at time t^{n+1} and psi_i at time t^* simulataneously: ------------------------#
    WF_h_psi = DW_solvers.WF_h_SE(dim, n_z, g, H, H0, Lw, WM, dWM_dy, dWM_dt, dt, delta_psi, delta_hat_star, h_n0, h_n1, x_coord, psi_1_n0, hat_psi_star, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN)
    
    #----------------------------------------- Step 2 : Update psi_1 at time t^{n+1}: -----------------------------------------#
    A_psi_s, L_psi_s = DW_solvers.WF_psi_SE(dim, g, H, H0, Lw, WM, WM_n1, dWM_dy, dWM_dt, dt, x_coord, delta_h, psi_1, psi_1_n0, hat_psi_star, h_n1, M11, MN1, MNN, D11, D1N, DN1, DNN,S11, SN1, SNN, A11, AN1, ANN, I1, IN)
    
    #----------------------------------------- Step 3 : Update psi_i at time t^{n+1}: -----------------------------------------#
    A_hat, L_hat = DW_solvers.WF_hat_psi_SE(dim, H, H0, g, n_z, Lw, x_coord, WM, dWM_dt, dWM_dy, dt, delta_hat_psi, hat_psi, h_n0, psi_1_n0, M11, MN1, MNN, D11, D1N, DN1, DNN,S11, SN1, SNN, A11, AN1, ANN, I1, IN)

elif scheme=="SV":#______________________________________________ Stormer-Verlet ______________________________________________#
    #--------------------------------------- Step 1 : Update psi_1^{n+1/2} and psi_i^*: ---------------------------------------#
    WF_psi_star = DW_solvers.WF_psi_half_SV(dim, n_z, g, H, H0, Lw, x_coord, WM, WM_half, dWM_dy, dWM_dt, dWM_half_dy, dWM_half_dt, dt, delta_psi, delta_hat_star, psi_1_n0, psi_1_half, hat_psi_star, h_n0, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN)
    #----------------------------- Step 2 : Update h^{n+1} and psi_i at time t^** simulataneously: ----------------------------#
    WF_h_psi = DW_solvers.WF_h_SV(dim, n_z, Lw, H0, g, dt, x_coord, WM, WM_half, dWM_half_dy, dWM_half_dt, delta_psi, delta_hat_star, h_n0, h_n1, psi_1_half, hat_psi_star, hat_psi_aux, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN)

    #----------------------------------------- Step 3 : Update psi_1 at time t^{n+1}: -----------------------------------------#
    a_psi_1, L_psi_1 = DW_solvers.WF_psi_n1_SV(dim, H0, H, g, x_coord, delta_h, Lw, WM_n1, WM_half, dt, psi_1_half, psi_1, dWM_half_dt, dWM_half_dy, hat_psi_aux, h_n1, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN )

    #----------------------------------------- Step 4 : Update psi_i at time t^{n+1}: -----------------------------------------#
    A_hat, L_hat = DW_solvers.WF_hat_psi_SV(dim, n_z, Lw, H0, H, WM, x_coord, dt, dWM_dt, dWM_dy, delta_hat_psi, hat_psi, h_n0, psi_1_n0, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN )

"""
    **************************************************************************************
    *                                 Define the solvers                                 *
    ************************************************************************************** """
#____________________________________ Solvers parameters ____________________________________#
param_h = {"ksp_converged_reason":False,"pc_type": "fieldsplit","pc_fieldsplit_type": "schur","pc_fieldsplit_schur_fact_type": "upper"}
param_psi = {"ksp_converged_reason":False,'ksp_type': 'preonly', 'pc_type': 'lu'}
param_hat_psi = {"ksp_converged_reason":False,'ksp_type': 'preonly', 'pc_type': 'lu'}

#--------------------------------------------------------------------------------------------#
#                                      Symplectic-Euler                                      #
#____________________________________________________________________________________________#
if scheme=="SE":
    #_______________________ Variational solver for h (and hat_psi^*) _______________________#
    h_problem = NonlinearVariationalProblem(WF_h_psi, w_n1)
    h_solver = NonlinearVariationalSolver(h_problem, solver_parameters=param_h)

    #_____________________________ Variational solver for psi_1 _____________________________#
    psi_problem = LinearVariationalProblem(A_psi_s, L_psi_s, psi_1_n1)
    psi_solver = LinearVariationalSolver(psi_problem, solver_parameters=param_psi)

    #____________________________ Variational solver for hat_psi ____________________________#
    hat_psi_problem = LinearVariationalProblem(A_hat, L_hat, hat_psi_n0)
    hat_psi_solver = LinearVariationalSolver(hat_psi_problem, solver_parameters=param_hat_psi)

#--------------------------------------------------------------------------------------------#
#                                       Stormer-Verlet                                       #
#____________________________________________________________________________________________#
if scheme=="SV":
    #_______________________ Variational solver for psi_1^{n+1/2} (and hat_psi^*) _______________________#
    psi_half_problem = NonlinearVariationalProblem(WF_psi_star, w_half)
    psi_half_solver = NonlinearVariationalSolver(psi_half_problem, solver_parameters=param_h)
    
    #_____________________________ Variational solver for h^{n+1} psi_i^** _____________________________#
    h_problem = NonlinearVariationalProblem(WF_h_psi, w_n1)
    h_solver = NonlinearVariationalSolver(h_problem, solver_parameters=param_h)
    
    #_______________________ Variational solver for psi_1^{n+1} _______________________#
    psi_n1_problem = LinearVariationalProblem(a_psi_1, L_psi_1, psi_1_n1)
    psi_n1_solver = LinearVariationalSolver(psi_n1_problem, solver_parameters=param_psi)
    
    #____________________________ Variational solver for hat_psi ____________________________#
    hat_psi_problem = LinearVariationalProblem(A_hat, L_hat, hat_psi_n0)
    hat_psi_solver = LinearVariationalSolver(hat_psi_problem, solver_parameters=param_hat_psi)

print('...solvers initialised')

"""
    *************************************************************
    *                        Saving Files                       *
    ************************************************************* """
save_waves, save_WM, WM_file, Energy_file, README_file = saving_files(save_path)

#_________ For measured data, create a file for each probe __________#
if input_data == "measurements":
    # Find the probe location and create the files:
    Ind_1, Ind_2, Ind_3, Ind_4, Ind_5, Ind_6, x1_file, x2_file, \
    x3_file, x4_file, x5_file,x6_file,all_file = probe_location(res, save_path)


"""
    ****************************************************************************
    *                                 Saving mesh                              *
    ****************************************************************************"""
#---------------------------------------------------------------------------------#
#                      Save waves in the 3D free-surface domain                   #
#---------------------------------------------------------------------------------#

if dim=='2D': # Extend the 1D horizontal mesh (x) to 2D horizontal mesh (x,y)
    mesh_2D = RectangleMesh(Nx,1,Lx,Ly,quadrilateral=True)        # 2D surface mesh
    V_2D = FunctionSpace(mesh_2D,"CG",1)                  # 2D surface funct. space
    Vec_2D = VectorFunctionSpace(mesh_2D,"CG",1, dim=n_z)  # 2D vector funct. space
    h_2D = Function(V_2D)                                                  # h(x,y)
    psi_s_2D = Function(V_2D)                                         # psi_1 (x,y)
    psi_i_2D = Function(Vec_2D)                                       # psi_i (x,y)
    beach_s_2D = Function(V_2D).interpolate(Expression("0.5*(1+copysign(1.0,x[0]-xb))*slope*(x[0]-xb)",
                                                       H0=H0,xb=xb, slope=sb))
    # b(x,y)
    # Extend the surface mesh in depth to obtain {0<x<Lx; 0<y<Ly; 0<z<H0}
    mesh_3D = ExtrudedMesh(mesh_2D,                   # horizontal mesh to extrude;
                           n_z,               # number of elements in the vertical;
                           layer_height=H0/(n_z),         # length of each element;
                           extrusion_type='uniform')     # type of extruded coord.;

else:# If the solutions are already (x,y)-dependent, we extend the domain in depth:
    mesh_3D = ExtrudedMesh(hor_mesh,                  # horizontal mesh to extrude;
                           n_z,               # number of elements in the vertical;
                           layer_height=H0/(n_z),         # length of each element;
                           extrusion_type='uniform')     # type of extruded coord.;


"""
    *****************************
    *      Function to save     *
    ***************************** """
#__________ Function Space _________#
V_3D = FunctionSpace(mesh_3D, "CG",1)
#____________ Functions ____________#
waves = Function(V_3D,name="phi")
WM_3D = Function(V_3D,name = "WM")

"""
    **************************************************************************
    *                         Mapping and transforms                         *
    **************************************************************************"""
if dim=="2D":
    # Indices to map h(x) and phi(x) to h(x,y) and phi(x,y) :
    Indx = []
    for j in range(len(hor_mesh.coordinates.dat.data[:])):
        Indx.append([y for y in range(len(mesh_2D.coordinates.dat.data[:,0]))\
         if mesh_2D.coordinates.dat.data[y,0]==hor_mesh.coordinates.dat.data[j]])

# Index used to differentiate each vertical layer
Indz = []
for i in range(0,n_z+1):
    Indz.append([zz for zz in range(len(mesh_3D.coordinates.dat.data[:,2])) \
     if mesh_3D.coordinates.dat.data[zz,2] == mesh_3D.coordinates.dat.data[i,2]])

# Index of the 3D funct. for which x<Lw. This is used to transform the 3D domain
# in x, to get back to the moving domain:
Test_x_Lw=Function(V_3D)
Test_x_Lw.interpolate(Expression('0.5*(1.0+copysign(1.0,Lw-x[0]))',Lw=Lw))
Indw = [item for item in range(len(Test_x_Lw.dat.data[:])) \
        if Test_x_Lw.dat.data[item] != 0.0]


print('Update of the solutions:')
""" *********************************************************************************
    *                                   Time loop                                   *
    ********************************************************************************* """
t_save = t
ind=0
it_time_save=0
run_time_it = time.time()-start_time
t_aux = t
update_wm = 'Yes' ;
while t<Tend-dt:
    """ *****************************************************************************
        *                               SAVE FUNCTIONS                              *
        ***************************************************************************** """
    if t_save <= t:
        print 'Progress: ', 100*t/Tend, ' %'
        #-------------------------------------------------------------------------------#
        #                                    ENERGY                                     #
        #-------------------------------------------------------------------------------#
        save_energy(h_n0, psi_1_n0, hat_psi_n0, WM, dWM_dy, dWM_dt, H, x_coord, Lw,
                    H0, g, A11, AN1, A1N, ANN, M11, M1N, MN1, MNN, D11, D1N, DN1,
                    DNN, S11, S1N, SN1, SNN, I1, IN, Energy_file, t, dim)
                    
                    
        #-------------------------------------------------------------------------------#
        #                               SAVE 3D FUNCTIONS                               #
        #-------------------------------------------------------------------------------#
        #______________________________ Project solutions ______________________________#
        if dim == '2D':
            # To the surface plane (x,y) :
            x_to_xy(h_n0, psi_1_n0, hat_psi_n0, h_2D, psi_s_2D, psi_i_2D, Indx)
            # In depth (x,y,z):
            for i in range(0,n_z+1):                                     # for each layer
                phi_projection(i, n_z, waves, Indz, psi_s_2D, psi_i_2D)  # phi(z) = psi_i
                WM_3D.dat.data[Indz[i]] = WM.dat.data[0]                     # WM(z) = WM
        elif dim == '3D':
            # In depth (x,y,z):
            for i in range(0,n_z+1):                                     # for each layer
                phi_projection(i, n_z, waves, Indz, psi_1_n0, hat_psi_n0)# phi(z) = psi_i
                WM_3D.dat.data[Indz[i]] = WM.dat.data[:]                     # WM(z) = WM

        #__________________________ Save the fixed coordinates _________________________#
        init_coord = mesh_3D.coordinates.vector().get_local()

        #_________________________________ z-transform _________________________________#
        if dim == '2D':
            z_transform(mesh_3D, n_z, h_2D, beach_s_2D, H0, Indz)
        elif dim == '3D':
            z_transform(mesh_3D, n_z, h_n0, beach, H0, Indz)

        #_________________________________ x-transform _________________________________#
        x_transform(mesh_3D, Lw, WM_3D, Indw)

        #_________________________________ Save waves __________________________________#
        save_waves.write(waves)

        #__________________________ Back to the initial mesh ___________________________#
        mesh_3D.coordinates.vector().set_local(init_coord)
        
        #_______________________________ Save wavemaker ________________________________#
        save_WM.write(WM_3D)
        

        #-------------------------------------------------------------------------------#
        #                        SURFACE DEVIATION AT THE PROBES                        #
        #-------------------------------------------------------------------------------#
        if input_data=="measurements":
            save_probes(t, h_n0, beach, Ind_1, Ind_2, Ind_3, Ind_4, Ind_5, Ind_6,
                        x1_file, x2_file, x3_file, x4_file, x5_file, x6_file,all_file)
            save_wm(t, WM, dWM_dt, WM_file, dWM_dt_file)

        #_____________________________ Update saving time ______________________________#
        t_save+=dt_save
        

    """ *********************************************************************
        *                            Update time                            *
        ********************************************************************* """

    #_______________________ Update time: t^n -> t^{n+1} _______________________#
    t_half = t+0.5*dt
    t += dt
    
    if input_data=="created":
        if t<=t_stop:                                    # The wavemaker keeps moving
            if scheme=="SV":
                WM_expr.t = t_half
                dWM_dt_expr.t = t_half
                dWM_dy_expr.t = t_half
                WM_half.interpolate(WM_expr)                              # update R(x,y;t)
                dWM_half_dt.interpolate(dWM_dt_expr)                         # update dR/dt
                dWM_half_dy.interpolate(dWM_dy_expr)                         # update dR/dy
            
            WM_expr.t = t
            dWM_dt_expr.t = t
            dWM_dy_expr.t = t
            WM_n1.interpolate(WM_expr)                              # update R(x,y;t)
            dWM_n1_dt.interpolate(dWM_dt_expr)                         # update dR/dt
            dWM_n1_dy.interpolate(dWM_dy_expr)                         # update dR/dy
            t_aux = t

        elif t>t_stop and update_wm=='Yes':# We stop the wavemaker motion;
            update_wm = 'No'
            if scheme=="SV":
                if t_half<=t_stop:
                    t_aux = t_half
            WM_expr.t = t_aux
            dWM_dt_expr.t = t_aux
            dWM_dy_expr.t = t_aux
            WM_n1.interpolate(WM_expr)
            dWM_n1_dt.assign(0.0)
            dWM_n1_dy.interpolate(dWM_dy_expr)
            if scheme=="SV":
                WM_half.interpolate(WM_expr)
                dWM_half_dt.assign(0.0)
                dWM_half_dy.interpolate(dWM_dy_expr)

    """ **************************************************
        *            Solve the weak formulations         *
        ************************************************** """
    #___________________ Call the solvers ___________________#
    if scheme=="SE":                     # 1st-order SE scheme
        h_solver.solve()           # get h^{n+1} and hat_psi^*
        psi_solver.solve()                     # get psi^{n+1}
    elif scheme=="SV":                   # 2nd-order SV scheme
        psi_half_solver.solve()# get psi^{n+1/2} and hat_psi^*   # WM_half
        h_solver.solve()        # get h^{n+1} and hat_psi^{**}   # WM_half
        psi_n1_solver.solve()                  # get psi^{n+1}   # WM_half
    """ *************************************************
        *               Update the functions            *
        ************************************************* """
    #_________________ Update the solutions ________________#
    h_out, hat_psi_out = w_n1.split()
    h_n0.assign(h_out)
    psi_1_n0.assign(psi_1_n1)
    hat_psi_n0.assign(hat_psi_out)

    #_________________ Update the wavemaker ________________#
    WM.assign(WM_n1)
    dWM_dt.assign(dWM_n1_dt)
    dWM_dy.assign(dWM_n1_dy)


comp_time = time.time()-start_time
jours = int(comp_time/(24*3600))
heures = int((comp_time-jours*24*3600)/3600)
minutes = int((comp_time-jours*24*3600-heures*3600)/60)
secondes = (comp_time -jours*24*3600-heures*3600 - minutes*60)/60
save_README(README_file, Lx, Ly, H0, xb, sb, res_x, Nx, Ny, Nz, gamma, Tw, w, t_stop, Lw, scheme, dt, t, jours, heures, minutes, secondes, comp_time)

