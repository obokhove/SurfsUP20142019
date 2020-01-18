"""
Coupling between deep- (potential flow) shallow- (nonlinear equations) models with Symplectic-Euler time scheme.
Author: Floriane Gidel. mmfg@leeds.ac.uk
"""

import os.path
import time
import numpy as np
from vertical_discr import *
import NLDW_WM as dw
import NLSW_beach as sw
from Settings import *
from measured_data import *
from savings import *
from firedrake import *



start_time = time.time()

"""
    *******************************************************************************
    *                                  Settings                                   *
    ******************************************************************************* """

#_____________________________________ Test case: ____________________________________#
input_data, test_case = input()
#_________________________________ Deep-water domain: ________________________________#
H0, xb, slope, H_expr, Hc, xc, Lw, res_dw, n_z = domain(input_data)

#_______________________________ Shallow-water domain: _______________________________#
Ldw = xc
Lsw = (1.2*H0 - (H0-Hc))/slope
res_sw = res_dw*res_dw

#___________________________________ Total domain ___________________________________#
L_total = Ldw + Lsw

#_______________________________________ Time: ______________________________________#
T0, t, dt, Tend, dt_save = set_time()

if input_data == "created":
    #__________________________________ Save path ___________________________________#
    save_path= 'numerical data/'+input_data+'/'+ test_case + '/'
    #__________________________________ Wavemaker: __________________________________#
    g, lamb, k, w, Tw, gamma, t_stop, WM_expr, dWM_expr = wavemaker(H0, Lw, t)

elif input_data == "measurements":
    g= 9.81
    #______________________ Run corresponding to the test_case ______________________#
    run_nr = run_number(test_case)
    save_path= 'numerical data/'+input_data+'/case '+test_case+'/'+run_nr+'/'
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    #_________________________________ Load the data ________________________________#
    measurement_path = 'measured data/case '+test_case + '/' + run_nr + '/'
    wm_data, wm_vel_data, t_data = load_wavemaker(measurement_path)

    #____________________________ Wavemaker expressions _____________________________#
    WM_expr, dWM_expr = interpolate_wavemaker(wm_data, wm_vel_data, t_data, t, dt, Lw)
    
    #_________________________________ Define time __________________________________#
    ind_data = 0
    T0 = t_data[ind_data]
    Tend = t_data[-1]
    t = T0
    t_stop = Tend
    
    #_________________ Indices corresponding to the probe locations _________________#
    Ind_1, Ind_2, Ind_3, Ind_4, Ind_5, Ind_6, Ind_7 = probe_location(res_dw)
    probe_files(save_path)

"""
    *******************************************************
    *                    Saving files                     *
    *******************************************************"""
dw_beach_file = File(os.path.join(save_path, "dw_beach.pvd"))
dw_h_save = File(os.path.join(save_path,"dw_h.pvd"))
sw_beach_file = File(os.path.join(save_path,"sw_beach.pvd"))
sw_h_file = File(os.path.join(save_path,"sw_h.pvd"))
dw_waves_file = File(os.path.join(save_path,"dw_waves.pvd"))
sw_waves_file = File(os.path.join(save_path,"sw_waves.pvd"))
WM_file = File(os.path.join(save_path,"wavemaker.pvd"))
E_file_dw = open(os.path.join(save_path,"energy_dw.txt"), 'w')
E_file_sw = open(os.path.join(save_path,"energy_sw.txt"), 'w')
README_file =open(os.path.join(save_path,"README.txt"), 'w')

if input_data == "measurements":
    x1_file, x2_file, x3_file, x4_file, x5_file, x6_file, x7_file = probe_files(save_path)

print('Settings uploaded')
print('Creation of the meshes...')
"""
    *********************************************
    *           Definition of the mesh          *
    ********************************************* """
#____________ Vertical discretization _____________#
Nz = n_z+1         # Number of points in one element
#___________ Horizontal discretization ____________#
Ne_dw = round(Ldw/res_dw) # Number of elements in x
Nv_sw = round(Lsw/res_sw)  # Number of volumes in x

#________________ Deep-water mesh _________________#
dw_mesh = IntervalMesh(Ne_dw,Ldw)         # DW mesh

#_______________ Shallow-water mesh _______________#
sw_mesh = IntervalMesh(Nv_sw,Lsw)         # SW mesh

"""
    *********************************************
    *     Definition of the function spaces     *
    ********************************************* """
#_________________ For h and psi_1 _________________#
dw_V = FunctionSpace(dw_mesh, "CG",1)
#___________________ For hat_psi ___________________#
dw_Vec = VectorFunctionSpace(dw_mesh,"CG",1,dim=n_z)
dw_W = dw_V*dw_Vec
#_________________ For sw h and hu _________________#
sw_V = FunctionSpace(sw_mesh,"DG", 0)

"""
    **************************************************
    *          Definition of the functions           *
    ************************************************** """
#____________________ Test functions ____________________#
v = TestFunction(dw_V)
q = TestFunction(dw_Vec)
w = TestFunction(dw_W)
p,r = split(w)

#____________________ Trial functions ___________________#
psi_s = TrialFunction(dw_V)
hat_psi = TrialFunction(dw_Vec)


#_______________________ Functions ______________________#
#---------------------- Deep water: ---------------------#
# Time t^n :
h_n0 = Function(dw_V)
psi_s_n0 = Function(dw_V)
hat_psi_n0 = Function(dw_Vec)
WM = Function(dw_V)
dWM = Function(dw_V)
x_coord = Function(dw_V).interpolate(Expression("x[0]"))
# Time t^n+1 :
w_n1 = Function(dw_W)
h_n1, hat_psi_aux = split(w_n1)
psi_s_n1 = Function(dw_V)
hat_psi_n1 = Function(dw_Vec)
WM_n1 = Function(dw_V)
dWM_n1 = Function(dw_V)
H = Function(dw_V)

#-------------------- Shallow water: --------------------#
h_fv = Function(sw_V, name="h_fv")
hu_fv = Function(sw_V, name = "hu_fv")

#------------------ BCs for deep-water: -----------------#
hu_fe = Function(dw_V)


"""
    ********************************************************************************
    *                     Save waves in the 2D free-surface domain                 *
    ********************************************************************************"""
#-------------------------------------------------------------------------------------#
#                        Save waves in the 3D free-surface domain                     #
#-------------------------------------------------------------------------------------#
#____________________________________ Deep water: ____________________________________#
dw_mesh_2D = RectangleMesh(Ne_dw,1,Ldw,1.0,quadrilateral=True)           # 2D surface mesh
dw_V_2D = FunctionSpace(dw_mesh_2D,"CG",1)                       # 2D surface funct. space
dw_Vec_2D = VectorFunctionSpace(dw_mesh_2D,"CG",1, dim=n_z)          # 2D vector funct. space
dw_h_2D = Function(dw_V_2D)                                                          # h(x,y)
dw_psi_s_2D = Function(dw_V_2D)                                                 # psi_1 (x,y)
dw_psi_i_2D = Function(dw_Vec_2D)                                               # psi_i (x,y)

dw_mesh_3D = ExtrudedMesh(dw_mesh_2D, n_z, layer_height=H0/(n_z),extrusion_type='uniform')
dw_V_3D = FunctionSpace(dw_mesh_3D,"CG",1)  # function space
dw_waves = Function(dw_V_3D, name="phi_dw")
dw_h = Function(dw_V, name="dw_h")

#___________________________________ Shallow water: __________________________________#
sw_mesh_2D = RectangleMesh(Nv_sw,1,Lsw,1.0,quadrilateral=True)           # 2D surface mesh
sw_mesh_3D = ExtrudedMesh(sw_mesh_2D, 1, layer_height = H0, extrusion_type='uniform')
sw_V_3D = FunctionSpace(sw_mesh_3D, "DG", 0)
sw_waves = Function(sw_V_3D, name = "sw_u")

sw_h = Function(sw_V, name="sw_h")
sw_u = Function(sw_V, name="sw_u")

#____________________________________ Wave-maker: ____________________________________#
mesh_WM = BoxMesh(2, 2, 2, 0.2, 1.0, 1.3*H0)
V_WM = FunctionSpace(mesh_WM,"CG",1)
WM_save = Function(V_WM)

#_______________________________________ Beach _______________________________________#
dw_beach = Function(dw_V)
sw_beach = Function(sw_V)
dw_beach_2D = Function(dw_V_2D).interpolate(Expression("0.5*(1+copysign(1.0,x[0]-xb))*slope*(x[0]-xb)",
                                                       H0=H0,xb=xb, slope=slope))

"""
    **************************************************************************
    *                         Mapping and transforms                         *
    **************************************************************************"""
#________________ Index used to extend from 1D to 2D ______________#
Indx_dw = []
for j in range(len(dw_mesh.coordinates.dat.data[:])):
    Indx_dw.append(np.where(dw_mesh_2D.coordinates.dat.data[:,0]==dw_mesh.coordinates.dat.data[j]))

#____________ Index used to differentiate each vertical layer (2D ->3D)_________#
#--------------------------------- Deep water: ---------------------------------#
coord = 2
Indz_dw = []
for i in range(0,n_z+1):
    Indz_dw.append([item for item \
                 in range(len(dw_mesh_3D.coordinates.dat.data[:,2])) \
                 if dw_mesh_3D.coordinates.dat.data[item,2] == dw_mesh_3D.coordinates.dat.data[i,2]])


#-------------------------------- Shallow water: -------------------------------#
Indz_sw = []
for i in range(len(sw_mesh.coordinates.dat.data[:])):
    Indz_sw.append(np.where(sw_mesh_3D.coordinates.dat.data[:,0]==sw_mesh.coordinates.dat.data[i]))


#_____________________________ Index  for which x<Lw ____________________________#
Test_ind_LW_DW = Function(dw_V_3D)
Test_ind_LW_DW.interpolate(Expression('0.5*(1.0+copysign(1.0,Lw-x[0]))',Lw=Lw))
Ind_Lw_DW = [item for item in range(len(Test_ind_LW_DW.dat.data[:])) if Test_ind_LW_DW.dat.data[item] != 0.0]

#---- x transform
for i in range(len(sw_mesh_3D.coordinates.dat.data[:,0])):
    sw_mesh_3D.coordinates.dat.data[i,0]+=xc
ii = np.where(sw_mesh_3D.coordinates.dat.data[:,0]==np.max(sw_mesh_3D.coordinates.dat.data[:,0]))

print('...meshes created')

"""
    ***********************************************************************************
    *                          Initialisation of the Functions                        *
    ***********************************************************************************"""
print('Initalisation of the functions...')
#----------------------------------------------------------------------------------------#
#                                       Topography                                       #
#----------------------------------------------------------------------------------------#
#___________________________________ Deep-water beach ___________________________________#
H.interpolate(H_expr)                                                 # Depth at rest H(x)
dw_beach.interpolate(H0-H)                                                    # Beach b(x)
# Save:
dw_beach_file.write(dw_beach)

#_________________________________ Shallow-water beach __________________________________#
beach_sw_expr = Expression("slope*x[0]",slope=slope)
sw_beach.interpolate(beach_sw_expr)
H0_sw = H0-sw_beach.dat.data[0]-dw_beach.dat.data[-1]

#----------------------------------------------------------------------------------------#
#                                       Wavemaker                                        #
#----------------------------------------------------------------------------------------#
WM.interpolate(WM_expr)
dWM.interpolate(dWM_expr)
WM_n1.interpolate(WM_expr)
dWM_n1.interpolate(dWM_expr)

#----------------------------------------------------------------------------------------#
#                                          Depth                                         #
#----------------------------------------------------------------------------------------#
#______________________________________ Deep water ______________________________________#
h_n0.assign(H)                                                         # h(x,y;t=0) = H(x)
w_n1.sub(0).assign(H)                                                            # h^{n+1}
#_____________________________________ Shallow water ____________________________________#
h_fv.assign(H0_sw-sw_beach) # faux car ca cest la valeur a xc mais dans lelement 0  h est une moyenne.
h_fv.dat.data[np.where(h_fv.vector().get_local()<0)]=0

#----------------------------------------------------------------------------------------#
#                                        Velocity                                        #
#----------------------------------------------------------------------------------------#
#______________________________________ Deep water ______________________________________#
psi_s_n0.assign(0.0)                                                 # \psi_1(x,y;t=0) = 0
for i in range(0,n_z):
    hat_psi_n0.dat.data[:,i] = 0.0                                               # psi_i^n
    w_n1.sub(1).dat.data[:,i] = 0.0                                            # psi_i^{*}
#_____________________________________ Shallow water ____________________________________#
hu_fv.assign(0.0)

print('...functions initialised')


"""
    *******************************************
    *  Compute vertical deep-water matrices   *
    ******************************************* """
print('Assembling z-matrices...')
#_____________ Initialization ____________#
Nz = n_z+1
A = np.eye(Nz,Nz)*0.0
M = np.eye(Nz,Nz)*0.0
D = np.eye(Nz,Nz)*0.0
S = np.eye(Nz,Nz)*0.0
Ik = np.eye(Nz,1)*0.0
Gk = np.eye(Nz,1)*0.0
    
    #______ Evaluation of each integral ______#
for i in range(0,Nz):
    for j in range(0,Nz):
        A[i,j]=A_ij(i,j,n_z,H0)
        M[i,j]=M_ij(i,j,n_z,H0)
        D[i,j]=D_ij(i,j,n_z,H0)
        S[i,j]=S_ij(i,j,n_z,H0)
    Ik[i] = I_i(i,n_z,H0)
    Gk[i] = G_i(i,n_z,H0)

#______________ Submatrices ______________#
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

G1 = Gk[0,0]
GN = as_tensor(Gk[1:,0])

Ii = []
for i in range(n_z):
    Ii.append(1.0)
I3 = as_tensor(Ii)

print('End of assembling')


"""
    **********************************************************************************************
    *                                 Define the deep-water solvers                              *
    ********************************************************************************************** """
print('Definition of the solvers...')
param_h = {"ksp_converged_reason":False,"pc_type": "fieldsplit",\
    "pc_fieldsplit_type": "schur","pc_fieldsplit_schur_fact_type": "upper"}
param_psi = {"ksp_converged_reason":False}
param_hat_psi = {"ksp_converged_reason":False,'ksp_type': 'preonly', 'pc_type': 'lu'}

#_______________________ Variational solver for h (and hat_psi^*) _______________________#
DW_VP_h = dw.VP_h(Lw, WM, p, h_n1, h_n0, dt, H0, psi_s_n0, x_coord, dWM, hat_psi_aux, r, w_n1, hu_fe, A11, A1N, AN1, ANN, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, S1N, SN1, SNN, I1, IN, I3, n_z)
DW_solver_h = NonlinearVariationalSolver(DW_VP_h, solver_parameters=param_h)
    
#_____________________________ Variational solver for psi_1 _____________________________#
DW_VP_psi_s = dw.VP_psi_s(Lw, H0, g, dt, WM, WM_n1, dWM, v, x_coord, psi_s_n0, h_n1, psi_s_n1, psi_s, H, hat_psi_aux, A11, A1N, AN1, ANN, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, S1N, SN1, SNN, I1, IN, I3, G1, GN, hu_fe)
DW_solver_psi_s = LinearVariationalSolver(DW_VP_psi_s, solver_parameters=param_psi)

#____________________________ Variational solver for hat_psi ____________________________#
DW_VP_hat_psi_n1 = dw.VP_hat_psi_n0(n_z, Lw, H0, g, dt, WM_n1, dWM_n1, q, h_n0, psi_s_n1,  hat_psi_n1, hat_psi, hu_fe,  A11, A1N, AN1, ANN, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, S1N, SN1, SNN, I1, IN, I3)
DW_solver_hat_psi_n1 = LinearVariationalSolver(DW_VP_hat_psi_n1, solver_parameters=param_hat_psi)

print('...solver defined')
print('Start of time iterations:')
"""
    ********************************************
    *        Solve the weak formulations       *
    ******************************************** """

t_save = t
E_sw=0
while t<Tend:
    print 'Progress: ' , 100*t/Tend, '%'
    """ *****************************************************************************
        *                               SAVE FUNCTIONS                              *
        ***************************************************************************** """
    if t_save<=t:
        print 'saving...'
        t_save+=dt_save
        #-------------------------------------------------------------------------------#
        #                                    ENERGY                                     #
        #-------------------------------------------------------------------------------#
        dw.energy(t, H, E_file_dw, Lw, H0, g, WM, dWM, h_n0, psi_s_n0, hat_psi_n0, A11, A1N, AN1, ANN, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, S1N, SN1, SNN, I1, IN, I3)
        E_file_sw.write('%-10s %-10s\n' % (str(t), str(E_sw)))
        
        #-------------------------------------------------------------------------------#
        #                                    DEPTH                                      #
        #-------------------------------------------------------------------------------#
        #__________________________________ Deep-Water _________________________________#
        dw_h.assign(h_n0 + dw_beach)
        dw_h_save.write(dw_h)
        #_________________________________ Shallow-Water _______________________________#
        # The mesh is shifted to start at xc: x = \check{x}+x_c.
        init_mesh = sw_mesh.coordinates.vector().get_local()
        sw_mesh.coordinates.dat.data[:]+=xc
        sw_h.assign(sw_beach.dat.data[0]+dw_beach.dat.data[-1] + sw_beach + h_fv)
        sw_h_file.write(sw_h)
        sw_mesh.coordinates.vector().set_local(init_mesh)


        #-------------------------------------------------------------------------------#
        #                                    WAVES                                      #
        #-------------------------------------------------------------------------------#
        #__________________________________ Deep-Water _________________________________#

        x_to_xy(h_n0, psi_s_n0, hat_psi_n0, dw_h_2D, dw_psi_s_2D, dw_psi_i_2D, Indx_dw)

        for i in range(0,n_z+1):                                         # for each layer
            #---- Projection of hat_psi and psi_1 #
            phi_projection(i, n_z, dw_waves, Indz_dw, dw_psi_s_2D, dw_psi_i_2D)
        
        init_DW_coord = dw_mesh_3D.coordinates.vector().get_local()
        #---- z transform
        z_transform(dw_mesh_3D, n_z, dw_h_2D, dw_beach_2D, H0, Indz_dw)
        #---- x transform
        dw_mesh_3D.coordinates.dat.data[Ind_Lw_DW[:],0]*=(Lw-WM.dat.data[0])/Lw
        dw_mesh_3D.coordinates.dat.data[Ind_Lw_DW[:],0]+=WM.dat.data[0]
        
        #---- save waves
        dw_waves_file.write(dw_waves)
        #---- transform back
        dw_mesh_3D.coordinates.vector().set_local(init_DW_coord)


        #_________________________________ Shallow-Water _______________________________#
        init_SW_coord = sw_mesh_3D.coordinates.vector().get_local()
        #---- z transform
        for i in range(len(h_fv.dat.data[:])):
            #---- Compute u
            sw_waves.dat.data[i]=hu_fv.dat.data[len(h_fv.dat.data[:])-1-i]/h_fv.dat.data[len(h_fv.dat.data[:])-1-i]
            sw_mesh_3D.coordinates.dat.data[Indz_sw[i],2]*=h_fv.dat.data[i]
            sw_mesh_3D.coordinates.dat.data[Indz_sw[i],2]+=(sw_beach.dat.data[i]+dw_beach.dat.data[-1]+sw_beach.dat.data[0])
        sw_mesh_3D.coordinates.dat.data[ii,2] = slope*(L_total-xb)
        sw_waves_file.write(sw_waves)
        #---- transform back
        sw_mesh_3D.coordinates.vector().set_local(init_SW_coord)#
        
        init_WM_coord = mesh_WM.coordinates.vector().get_local()
        mesh_WM.coordinates.dat.data[:,0]+=WM.dat.data[0]-np.max(mesh_WM.coordinates.dat.data[:,0])
        WM_file.write(WM_save)
        mesh_WM.coordinates.vector().set_local(init_WM_coord)


        if input_data == "measurements":
            save_probes(t, h_n0, dw_beach, x1_file, x2_file, x3_file, x4_file, x5_file, x6_file, x7_file, Ind_1, Ind_2, Ind_3, Ind_4, Ind_5, Ind_6, Ind_7)

    """ *********************************************************************
        *                            Update time                            *
        ********************************************************************* """
    t+=dt

    if input_data=="measurements":
        ind_data, dWM_expr, WM_expr = update_input_data(ind_data, t, t_data, wm_data, \
                                                        wm_vel_data, dWM_expr, \
                                                        WM_expr)
    if t<=t_stop:
        WM_expr.t = t
        dWM_expr.t = t
        WM_n1.interpolate(WM_expr)
        dWM_n1.interpolate(dWM_expr)
    else:
        dWM_expr.t = t_stop
        dWM_n1.interpolate(dWM_expr)
        WM_n1.assign(WM)

    """ **************************************************
        *            Solve the weak formulations         *
        ************************************************** """
    #______________ Call the deep-water solvers ______________#
    DW_solver_h.solve()                      # h^{n+1}, psi_i^*
    DW_solver_psi_s.solve()                       # psi_1^{n+1}

    #________ Update the boundary conditions solutions _______#
    h_out, hat_psi_out = w_n1.split()

    hu_bc = assemble((1/H0)*(h_out*psi_s_n0.dx(0)*I1 \
                             + h_out*dot(hat_psi_out.dx(0),IN)\
                             - G1*psi_s_n0*h_out.dx(0) \
                             - h_out.dx(0)*dot(GN,hat_psi_out))*ds(2))
    h_bc = assemble((h_out)*ds(2))

    #______________ Call the shallow-water solver _____________#
    h_fv, hu_fv, U, hu_fe, E_sw = sw.solve_FV(int(Nv_sw), res_sw, dt, sw_beach, g, h_bc, hu_bc, h_fv, hu_fv, hu_fe)


    #___________ Update the solutions for next step ___________#
    h_n0.assign(h_out)
    psi_s_n0.assign(psi_s_n1)
    hat_psi_n0.assign(hat_psi_out)
    WM.assign(WM_n1)
    dWM.assign(dWM_n1)

comp_time = time.time()-start_time

print('End of simulations')
print 'Computational time: ', comp_time
if input_data =="measurements":
    save_README(README_file, input_data, Ldw, Lsw, L_total, H0, xb, slope, xc, Lw, res_dw, Ne_dw, res_sw, Nv_sw, Nz, 0.0, 0.0, 0.0, t_stop, test_case, run_nr, t, dt, comp_time)
elif input_data =="created":
    save_README(README_file, input_data, Ldw, Lsw, L_total, H0, xb, slope, xc, Lw, res_dw, Ne_dw, res_sw, Nv_sw, Nz, gamma, Tw, w, t_stop, 'none', 'none', t, dt, comp_time)







