from firedrake import *
import os.path

"""
    *************************************************************
    *                        Saving Files                       *
    ************************************************************* """
def saving_files(save_path):
    save_waves = File(os.path.join(save_path, "waves.pvd"))
    save_WM = File(os.path.join(save_path, "Wavemaker.pvd"))
    WM_file = open(os.path.join(save_path, 'wm_motion.txt'), 'w')
    Energy_file = open(os.path.join(save_path, 'energy.txt'), 'w')
    README_file = open(os.path.join(save_path, 'README.txt'), 'w')
    return save_waves, save_WM, WM_file, Energy_file, README_file



"""
    ***************************************************************************
    *                                 Saving mesh                             *
    ***************************************************************************"""
def saving_domain(dim, hor_mesh, Nx, Ny, n_z, Lx, Ly, H0):

    if dim=='2D': # Extend the 1D horizontal mesh (x) to 2D horizontal mesh (x,y)
        mesh_2D = RectangleMesh(Nx,1,Lx,Ly,quadrilateral=True)       # 2D surface mesh
        V_2D = FunctionSpace(mesh_2D,"CG",1)                 # 2D surface funct. space
        Vec_2D = VectorFunctionSpace(mesh_2D,"CG",1, dim=n_z) # 2D vector funct. space
        h_2D = Function(V_2D)                                                 # h(x,y)
        psi_s_2D = Function(V_2D)                                        # psi_1 (x,y)
        psi_i_2D = Function(Vec_2D)                                      # psi_i (x,y)
        beach_s_2D = Function(V_2D).interpolate(beach_expr)                   # b(x,y)
        # Extend the surface mesh in depth to obtain {0<x<Lx; 0<y<Ly; 0<z<H0}
        mesh_3D = ExtrudedMesh(mesh_2D,                  # horizontal mesh to extrude;
                               n_z,              # number of elements in the vertical;
                               layer_height=H0/(n_z),        # length of each element;
                               extrusion_type='uniform')    # type of extruded coord.;
            
        # Indices to map h(x), phi(x) to h(x,y) and phi(x,y) :
        Indx = []
        for j in range(len(hor_mesh.coordinates.dat.data[:])):
            Indx.append([item for item in range(len(mesh_2D.coordinates.dat.data[:,0]))\
            if mesh_2D.coordinates.dat.data[item,0]==hor_mesh.coordinates.dat.data[j]])

    else:#If the solutions are already (x,y)-dependent, we extend the domain in depth:
        mesh_3D = ExtrudedMesh(hor_mesh,                 # horizontal mesh to extrude;
                               n_z,              # number of elements in the vertical;
                               layer_height=H0/(n_z),        # length of each element;
                               extrusion_type='uniform')    # type of extruded coord.;



""" *******************************************************
    *                  Compute the energy                 *
    ******************************************************* """
def save_energy(h_n0, psi_1_n0, hat_psi_n0, WM, dWM_dy, dWM_dt, H, x_coord, Lw, H0, g, A11, AN1, A1N, ANN, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, S1N, SN1, SNN, I1, IN, Energy_file, t, dim):
    if dim =="3D":
        energy = assemble(( 0.5*(h_n0/(Lw-WM))*(Lw**2 \
                                                + ((x_coord-Lw)*dWM_dy)**2)*((psi_1_n0.dx(0)**2)*M11 \
                                                                             + dot(hat_psi_n0.dx(0), \
                                                                                   (2.0*MN1*psi_1_n0.dx(0) \
                                                                                    + dot(MNN,hat_psi_n0.dx(0))))) \
                           +0.5*(Lw-WM)*h_n0*((psi_1_n0.dx(1)**2)*M11 \
                                              + dot(hat_psi_n0.dx(1), \
                                                    (2.0*MN1*psi_1_n0.dx(1) \
                                                     + dot(MNN,hat_psi_n0.dx(1)))))\
                           +(x_coord-Lw)*dWM_dy*h_n0*(psi_1_n0.dx(0)*(M11*psi_1_n0.dx(1) + dot(M1N,hat_psi_n0.dx(1))) \
                                                      +dot(hat_psi_n0.dx(0), (MN1*psi_1_n0.dx(1) \
                                                                              +dot(MNN,hat_psi_n0.dx(1)))))\
                           -((1/(Lw-WM))*(Lw**2 + ((x_coord-Lw)*dWM_dy)**2)*h_n0.dx(0)\
                             + (x_coord-Lw)*dWM_dy*h_n0.dx(1))*(psi_1_n0.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_n0)) \
                                                                +dot(hat_psi_n0.dx(0), (DN1*psi_1_n0 \
                                                                                        +dot(DNN, hat_psi_n0))))\
                           -((Lw-WM)*h_n0.dx(1) + (x_coord-Lw)*dWM_dy*h_n0.dx(0))*(psi_1_n0.dx(1)*(D11*psi_1_n0 \
                                                                                                   +dot(D1N,hat_psi_n0)) \
                                                                                   +dot(hat_psi_n0.dx(1), \
                                                                                        (DN1*psi_1_n0 \
                                                                                         + dot(DNN, hat_psi_n0)))) \
                           +(1/h_n0)*((0.5/(Lw-WM))*(h_n0.dx(0)**2)*(Lw**2 + ((x_coord-Lw)*dWM_dy)**2)\
                                      +0.5*(Lw-WM)*(h_n0.dx(1)**2) \
                                      + (x_coord-Lw)*dWM_dy*h_n0.dx(0)*h_n0.dx(1))*(S11*psi_1_n0**2 \
                                                                                    + dot(hat_psi_n0, 2.0*SN1*psi_1_n0 \
                                                                                          + dot(SNN,hat_psi_n0)))\
                           +0.5*((Lw-WM)*H0*H0/h_n0)*(A11*psi_1_n0**2 + dot(hat_psi_n0, (2.0*AN1*psi_1_n0 \
                                                                                         + dot(ANN,hat_psi_n0))))\
                           +0.5*H0*g*(Lw-WM)*(h_n0-H)**2)*dx)
    elif dim=="2D":
        energy = assemble(( 0.5*(h_n0/(Lw-WM))*(Lw**2)*((psi_1_n0.dx(0)**2)*M11 \
                                                        + dot(hat_psi_n0.dx(0), \
                                                              (2.0*MN1*psi_1_n0.dx(0) \
                                                               + dot(MNN,hat_psi_n0.dx(0))))) \
                           -((1/(Lw-WM))*(Lw**2)*h_n0.dx(0))*(psi_1_n0.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_n0)) \
                                                              +dot(hat_psi_n0.dx(0), \
                                                                   (DN1*psi_1_n0+dot(DNN, hat_psi_n0))))\
                           +(1/h_n0)*((0.5/(Lw-WM))*(h_n0.dx(0)**2)*(Lw**2))*(S11*psi_1_n0**2 \
                                                                              + dot(hat_psi_n0, 2.0*SN1*psi_1_n0 \
                                                                                    + dot(SNN,hat_psi_n0)))\
                           +0.5*((Lw-WM)*H0*H0/h_n0)*(A11*psi_1_n0**2 + dot(hat_psi_n0, \
                                                                            (2.0*AN1*psi_1_n0 + dot(ANN,hat_psi_n0))))\
                           +0.5*H0*g*(Lw-WM)*(h_n0-H)**2)*dx)


    Energy_file.write('%-10s %-10s %-10s %-10s %-10s %-10s %-10s\n' % (str(t),'   ', str(energy), '    ' , str(h_n0.dat.data[0]), '   ', str(dWM_dt.dat.data[0])))



#----------------------------------------------------------------------#
#                        Surface solutions (x,y)                       #
#----------------------------------------------------------------------#
def x_to_xy(h_n0, psi_1_n0, hat_psi_n0, h_2D, psi_s_2D, psi_i_2D, Indx):
    for i in range(len(h_n0.dat.data[:])):
        h_2D.dat.data[Indx[i]] = h_n0.dat.data[i]
        psi_s_2D.dat.data[Indx[i]]=psi_1_n0.dat.data[i]
        psi_i_2D.dat.data[Indx[i],:] = hat_psi_n0.dat.data[i,:]

#------------------------------------------------------------------------#
#                           3D solution (x,y,z)                          #
#------------------------------------------------------------------------#
def phi_projection(i, n_z, waves, Indz, psi_s, psi_i):
    if i==n_z:                                                   # if i=1,
        waves.dat.data[Indz[i]] = psi_s.dat.data[:]       # phi(z_i)=psi_1
    else:                                                        # if i>1,
        waves.dat.data[Indz[i]] = psi_i.dat.data[:,n_z-1-i] # phi(z_i)=psi_i


#-------------------------------------------------------------------------------------#
#                                 Transform the domain                                #
#-------------------------------------------------------------------------------------#

#____________________________________ z-transform ____________________________________#
def z_transform(mesh_3D, n_z, h_2D, beach_2D, H0, Indz):
    for i in range(0, n_z+1):                                          # for each layer
        mesh_3D.coordinates.dat.data[Indz[i],2]*=h_2D.dat.data[:]/H0  # z -> z*h/H0
        mesh_3D.coordinates.dat.data[Indz[i],2]+=beach_2D.dat.data[:] # z -> z+b(x)

#____________________________________ x-transform ____________________________________#
def x_transform(mesh_3D, Lw, WM_3D, Indw):
    for i in range(0,len(Indw)): # x -> R + x*(Lw-R)/Lw
        mesh_3D.coordinates.dat.data[Indw[i],0]*=(Lw-WM_3D.dat.data[Indw[i]])/Lw
        mesh_3D.coordinates.dat.data[Indw[i],0]+=WM_3D.dat.data[Indw[i]]












def save_README(README_file, Lx, Ly, H0, xb, sb, res, Nx, Ny, Nz, gamma, Tw, w, t_stop, Lw, scheme, dt, t, jours, heures, minutes, secondes, comp_time):
    README_file.write('______________________________________\n')
    README_file.write('                Summary               \n')
    README_file.write('--------------------------------------\n\n')

    README_file.write('------ Dimensions of the domain ------\n')
    README_file.write('Lenght Lx: %-5s m\n' %(str(Lx)))
    README_file.write('Lenght Ly: %-5s m\n' %(str(Ly)))
    README_file.write('Depth H0: %-5s m\n' %(str(H0)))
    README_file.write('Beach start: %-5s m\n' %str(xb))
    README_file.write('Beach slope: %-5s\n\n' %str(sb))

    README_file.write('----------- Mesh resolution ----------\n')
    README_file.write('In x: %-5s m (%-5s elements)\n' % (str(res),str(Nx)))
    README_file.write('In y: %-5s m (%-5s elements)\n' % (str(res),str(Ny)))
    README_file.write('In z: %-5s m (%-5s elements)\n\n' % (str(H0/Nz),str(Nz)))

    README_file.write('-------------- Wavemaker -------------\n')
    README_file.write('Amplitude: %-5s m\n' %(str(gamma)))
    README_file.write('Period: %-5s s\n' %(str(Tw)))
    README_file.write('Frequency: %-5s /s\n' %(str(w)))
    README_file.write('Stops after %-5s periods ( = %-5s s)\n' %(str(t_stop/Tw), str(t_stop)))
    README_file.write('Lw: %-5s m\n\n' %(str(Lw)))

    README_file.write('--------------- Solver ---------------\n')
    if scheme=="SE":
        README_file.write('1st order Symplectic-Euler scheme\n\n')
    else:
        README_file.write('2nd order Stormer-Verlet scheme\n\n')

    README_file.write('------------- Final time -------------\n')
    README_file.write('Tend: %-5s\n' %(str(t)))
    README_file.write('dt: %-5s\n' %(str(dt)))
    README_file.write('Computational time: %-5s j %-5s h %-5s mn %-5s s (=%-10s s)\n' %(str(jours),str(heures),str(minutes), str(secondes),str(comp_time)))


""" **************************************
    *   Updated interpolation with time  *
    **************************************"""
def update_input_data(ind, t, t_data, wm_data, \
                      wm_vel_data, dWM_dt_expr, \
                      WM_expr):
    #  if t<t2, update the interpolated interval:
    if t_data[ind+1]<=t:           # if t>t2:
        ind +=1
        #--------- update t1 and t2 --------#
        dWM_dt_expr.t1 = t_data[ind]
        dWM_dt_expr.t2 = t_data[ind+1]
        WM_expr.t1 = t_data[ind]
        WM_expr.t2 = t_data[ind+1]
        #----- update the velocity data ----#
        dWM_dt_expr.dwm1 = wm_vel_data[ind]
        dWM_dt_expr.dwm2 = wm_vel_data[ind+1]
        #------ update the motion data -----#
        WM_expr.wm1 = wm_data[ind]
        WM_expr.wm2 = wm_data[ind+1]

    return ind, dWM_dt_expr, WM_expr

""" ***********************************************
    *               Probe locations               *
    ***********************************************"""
def probe_location(res, save_path):
    
    x1 = 10
    x2 = 20
    x3 = 40
    x4 = 49.5
    x5 = 50
    x6 = 54
    
    #- Indices corresponding to the probe locations -#
    Ind_1 = int(x1/res)
    Ind_2 = int(x2/res)
    Ind_3 = int(x3/res)
    Ind_4 = int(x4/res)
    Ind_5 = int(x5/res)
    Ind_6 = int(x6/res)
    
    #------- Create saving file for each probe ------#
    x1_file = open(os.path.join(save_path, 'probe1.txt'), 'w')
    x2_file = open(os.path.join(save_path, 'probe2.txt'), 'w')
    x3_file = open(os.path.join(save_path, 'probe3.txt'), 'w')
    x4_file = open(os.path.join(save_path, 'probe4.txt'), 'w')
    x5_file = open(os.path.join(save_path, 'probe5.txt'), 'w')
    x6_file = open(os.path.join(save_path, 'probe6.txt'), 'w')
    all_file = open(os.path.join(save_path, 'all_probes.txt'), 'w')

    return Ind_1, Ind_2, Ind_3, Ind_4, Ind_5, Ind_6, \
        x1_file, x2_file, x3_file, x4_file, x5_file, x6_file, all_file

""" ******************************************************************
    *    Save the free-surface elevation at the probes' locations    *
    ****************************************************************** """
def save_probes(t, h, beach, Ind_1, Ind_2, Ind_3, Ind_4, Ind_5, Ind_6, \
                   x1_file, x2_file, x3_file, x4_file, x5_file, x6_file, all_file):
    

    #---------------- wave elevation at all probes -----------------#
    all_file.write('%-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n'
                  %(str(t), '   ',
                    str(h.dat.data[Ind_1]-beach.dat.data[Ind_1]), '   ',
                    str(h.dat.data[Ind_2]-beach.dat.data[Ind_2]), '   ',
                    str(h.dat.data[Ind_3]-beach.dat.data[Ind_3]), '   ',
                    str(h.dat.data[Ind_4]-beach.dat.data[Ind_4]), '   ',
                    str(h.dat.data[Ind_5]-beach.dat.data[Ind_5]), '   ',
                    str(h.dat.data[Ind_6]-beach.dat.data[Ind_6])))
    #---------------- wave elevation probe 1 : x1 = 10m -----------------#
    x1_file.write('%-10s \n' %(str(h.dat.data[Ind_1]-beach.dat.data[Ind_1])))
    #---------------- wave elevation probe 2 : x2 = 20m -----------------#
    x2_file.write('%-10s \n' %( str(h.dat.data[Ind_2]-beach.dat.data[Ind_2])))
    #---------------- wave elevation probe 3 : x3 = 40m -----------------#
    x3_file.write('%-10s\n' %( str(h.dat.data[Ind_3]-beach.dat.data[Ind_3])))
    #--------------- wave elevation probe 4 : x4 = 49.5m ----------------#
    x4_file.write('%-10s \n' %( str(h.dat.data[Ind_4]-beach.dat.data[Ind_4])))
    #---------------- wave elevation probe 4 : x4 = 50m -----------------#
    x5_file.write('%-10s \n' %(str(h.dat.data[Ind_5]-beach.dat.data[Ind_5])))
    #---------------- wave elevation probe 4 : x4 = 54m -----------------#
    x6_file.write('%-10s \n' %(str(h.dat.data[Ind_6]-beach.dat.data[Ind_6])))

""" ******************************************************************
    *             Save the wavemaker motion and velocity             *
    ****************************************************************** """
def save_wm(t, WM, dWM_dt, WM_file, dWM_dt_file):
    # Wave maker motion
    WM_file.write('%-10s %-10s\n' %(str(t), str(WM.dat.data[0])))
    # Wave maker velocity
    dWM_dt_file.write('%-10s %-10s\n' %(str(t), str(dWM_dt.dat.data[0])))


