from firedrake import *
import os.path



#-------------------------------------------------------------------------------------#
#                                  3D solution (x,y,z)                                #
#-------------------------------------------------------------------------------------#
#_____________________________________ projection ____________________________________#
#----------------------------------------------------------------------#
#                        Surface solutions (x,y)                       #
#----------------------------------------------------------------------#
def x_to_xy(h_n0, psi_1_n0, hat_psi_n0, h_2D, psi_s_2D, psi_i_2D, Indx):
    for i in range(len(h_n0.dat.data[:])):
        h_2D.dat.data[Indx[i]] = h_n0.dat.data[i]
        psi_s_2D.dat.data[Indx[i]]=psi_1_n0.dat.data[i]
        psi_i_2D.dat.data[Indx[i],:] = hat_psi_n0.dat.data[i,:]

def phi_projection(i, n_z, waves, Indz, psi_s, psi_i):
    if i==n_z:                                                   # if i=1,
        waves.dat.data[Indz[i]] = psi_s.dat.data[:]       # phi(z_i)=psi_1
    else:                                                        # if i>1,
        waves.dat.data[Indz[i]] = psi_i.dat.data[:,n_z-1-i] # phi(z_i)=psi_i

#____________________________________ z-transform ____________________________________#
def z_transform(mesh_3D, n_z, h_2D, beach_2D, H0, Indz):
    for i in range(0, n_z+1):                                          # for each layer
        mesh_3D.coordinates.dat.data[Indz[i],2]*=h_2D.dat.data[:]/H0  # z -> z*h/H0
        mesh_3D.coordinates.dat.data[Indz[i],2]+=beach_2D.dat.data[:] # z -> z+b(x)

#------------------------------------------------------------------------#
#                               README FILE                              #
#------------------------------------------------------------------------#

def save_README(README_file, input_data, Ldw, Lsw, L, H0, xb, slope, xc, Lw, res_dw, Ne_dw, res_sw, Nv_sw,  Nz, gamma, Tw, w, t_stop, test_case, run_nr, t, dt, comp_time):
    README_file.write('______________________________________\n')
    README_file.write('                Summary               \n')
    README_file.write('--------------------------------------\n\n')

    README_file.write('------ Dimensions of the domain ------\n')
    README_file.write('Deep-water lenght Lx: %-5s m\n' %(str(Ldw)))
    README_file.write('Shallow-water lenght Ly: %-5s m\n' %(str(Lsw)))
    README_file.write('Total lenght Ly: %-5s m\n' %(str(L)))
    README_file.write('Depth H0: %-5s m\n' %(str(H0)))
    README_file.write('Beach start: %-5s m\n' %str(xb))
    README_file.write('Beach slope: %-5s\n\n' %str(slope))
    README_file.write('Coupling at xc = %-5s m\n' %(str(xc)))
    README_file.write('x-transfrom on %-5s m\n' %(str(Lw)))
    
    README_file.write('----------- Mesh resolution ----------\n')
    README_file.write('In deep water: %-5s m (%-5s elements)\n' % (str(res_dw),str(Ne_dw)))
    README_file.write('In shallow water: %-5s m (%-5s elements)\n' % (str(res_sw),str(Nv_sw)))
    README_file.write('In z: %-5s m (%-5s elements)\n\n' % (str(H0/Nz),str(Nz)))
    
    if input_data =="created":
        README_file.write('-------------- Wavemaker -------------\n')
        README_file.write('Amplitude: %-5s m\n' %(str(gamma)))
        README_file.write('Period: %-5s s\n' %(str(Tw)))
        README_file.write('Frequency: %-5s /s\n' %(str(w)))
        README_file.write('Stops after %-5s periods ( = %-5s s)\n' %(str(t_stop/Tw), str(t_stop)))
    elif input_data == "measurements":
        README_file.write('-------------- Wavemaker -------------\n')
        README_file.write('Wavemaker input imported from measured data, in test case %-5s, run number %-5s' %(test_case, run_nr))

    README_file.write('--------------- Solver ---------------\n')
    README_file.write('1st order Symplectic-Euler scheme\n\n')

    README_file.write('------------- Final time -------------\n')
    README_file.write('Tend: %-5s\n' %(str(t)))
    README_file.write('dt: %-5s\n' %(str(dt)))
    README_file.write('Computational time: %-5s s \n' %(str(comp_time)))


