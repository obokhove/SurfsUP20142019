"""
    Solving 2D nonlinear potential flow equations with wave maker"
    Created 26 May 2016"
    Author : Floriane Gidel"
    """
import time
from vertical_discr import *
from firedrake import *


"""
    *********************************************
    *           Definition of the mesh          *
    ********************************************* """

"""
    ********************************************
    *       Define the weak formulations       *
    ******************************************** """

def VP_hat_psi_n0(n_z, Lw, H0, g, dt, WM, dWM_dt, delta_hat_psi, h_n0, psi_1_n0, hat_psi_n0, hat_psi, hu_fe,  A11, A1N, AN1, ANN, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, S1N, SN1, SNN, I1, IN, I3):
    a_hat_psi =((h_n0/(Lw-WM))*(Lw*Lw)*elem_mult(delta_hat_psi.dx(0), dot(MNN,hat_psi.dx(0)))\
                -((Lw*Lw)*h_n0.dx(0)/(Lw-WM))*(elem_mult(delta_hat_psi, dot(DNN.T,hat_psi.dx(0))) \
                                               + elem_mult(delta_hat_psi.dx(0),dot(DNN,hat_psi))) \
                +(1.0/h_n0)*((Lw*Lw)*(h_n0.dx(0)**2)/(Lw-WM))*elem_mult(delta_hat_psi,dot(SNN,hat_psi))\
                + ((Lw-WM)*H0*H0/h_n0)*elem_mult(delta_hat_psi,dot(ANN,hat_psi)))
        
        
    L_hat_psi =-((h_n0/(Lw-WM))*(Lw*Lw)*elem_mult(delta_hat_psi.dx(0), MN1*psi_1_n0.dx(0))\
                 -((Lw*Lw)*h_n0.dx(0)/(Lw-WM))*(elem_mult(delta_hat_psi, D1N*psi_1_n0.dx(0)) \
                                                + elem_mult(delta_hat_psi.dx(0),DN1*psi_1_n0)) \
                 +(1.0/h_n0)*( (Lw*Lw)*(h_n0.dx(0)**2)/(Lw-WM))*elem_mult(delta_hat_psi,SN1*psi_1_n0)\
                 + ((Lw-WM)*H0*H0/h_n0)*elem_mult(delta_hat_psi,AN1*psi_1_n0))
                 
                 
#    a_BC = -Lw*H0*elem_mult(delta_hat_psi, dot(DNN,hat_psi))
    hat_psi_BC = -(Lw*dWM_dt*h_n0*elem_mult(delta_hat_psi,IN))
    WF_hat_bound_2 = Lw*(hu_fe*elem_mult(delta_hat_psi,IN)) #+h_n0.dx(0)*(psi_1_n0*elem_mult(delta_hat_psi,DN1)))
    # PEUT ETRE ENLEVER LE H0
    
    A_hat = sum((a_hat_psi[ind])*dx for ind in range(0,n_z)) #+ sum((a_BC[ind])*ds(2) for ind in range(0,n_z))
    L_hat = sum((L_hat_psi[ind])*dx for ind in range(0,n_z)) + sum((hat_psi_BC[ind])*ds(1) for ind in range(0,n_z))+ sum((WF_hat_bound_2[ind])*ds(2) for ind in range(0,n_z))
    VP = LinearVariationalProblem(A_hat, L_hat, hat_psi_n0)

    return VP


def VP_h(Lw, WM, delta_psi, h_n1, h_n0, dt, H0, psi_1_n0, x_coord, dWM_dt, hat_psi_star, delta_hat_star, w1, hu_fe,  A11, A1N, AN1, ANN, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, S1N, SN1, SNN, I1, IN, I3, n_z):
    # Get h at time n+1. Implicit step for h_n1 : Nonlinear solver.

    # CHECK THE COUPLING CONDITIONS
    WF_h = (H0*delta_psi*(h_n1-h_n0)*(Lw-WM)/dt \
        -((h_n1/(Lw-WM))*(Lw*Lw)*(psi_1_n0.dx(0)*M11 + dot(hat_psi_star.dx(0),MN1))*delta_psi.dx(0)\
          -( (1/(Lw-WM))*(Lw*Lw)*h_n1.dx(0))*( delta_psi.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_star)) \
                                              +delta_psi*(psi_1_n0.dx(0)*D11 + dot(hat_psi_star.dx(0),DN1)))\
          +(1/h_n1)*((1/(Lw-WM))*(Lw*Lw)*(h_n1.dx(0)**2))*(psi_1_n0*S11 + dot(hat_psi_star,SN1))*delta_psi\
          +((Lw-WM)*H0*H0/h_n1)*(psi_1_n0*A11 + dot(hat_psi_star,AN1))*delta_psi \
          -delta_psi*H0*(x_coord-Lw)*dWM_dt*h_n1.dx(0)))*dx - (delta_psi*Lw*dWM_dt*h_n1*I1)*ds(1) + (delta_psi*(hu_fe*I1*Lw))*ds(2)
          
          
    WF_hat_psi_star= -((h_n1/(Lw-WM))*(Lw*Lw)*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_n0.dx(0)+dot(MNN,hat_psi_star.dx(0))))\
                      -((Lw*Lw)*h_n1.dx(0)/(Lw-WM))*(elem_mult(delta_hat_star, (psi_1_n0.dx(0)*D1N+ dot(DNN.T,hat_psi_star.dx(0)))) \
                                                     + elem_mult(delta_hat_star.dx(0),(DN1*psi_1_n0+dot(DNN,hat_psi_star)))) \
                      +(1.0/h_n1)*((Lw*Lw)*(h_n1.dx(0)**2)/(Lw-WM))*elem_mult(delta_hat_star,(SN1*psi_1_n0+ dot(SNN,hat_psi_star)))\
                      + ((Lw-WM)*H0*H0/h_n1)*elem_mult(delta_hat_star,(AN1*psi_1_n0+dot(ANN,hat_psi_star))))
        
        
    WF_hat_bound = -dWM_dt*elem_mult(IN,delta_hat_star)*h_n1*Lw
    WF_hat_bound_2 = Lw*(hu_fe*elem_mult(delta_hat_star,IN))#-h_n1.dx(0)*(psi_1_n0*elem_mult(delta_hat_star,DN1)+elem_mult(delta_hat_star, dot(DNN,hat_psi_star))))

    WF1 = WF_h + sum((WF_hat_psi_star[ind])*dx for ind in range(0,n_z)) + sum((WF_hat_bound[ind])*ds(1) for ind in range(0,n_z)) + sum((WF_hat_bound_2[ind])*ds(2) for ind in range(0,n_z))
    VP = NonlinearVariationalProblem(WF1, w1)
    return VP

###################

def VP_psi_s(Lw, H0, g, dt, WM, WM_n1, dWM_dt, delta_h, x_coord, psi_1_n0, h_n1, psi_1_n1, psi_1, H, hat_psi_star, A11, A1N, AN1, ANN, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, S1N, SN1, SNN, I1, IN, I3, G1, GN, hu_fe):

    A_psi_s = (H0*delta_h*(Lw-WM_n1)*psi_1)*dx
    
    L_psi_s = -(-H0*delta_h*(Lw-WM)*psi_1_n0 \
                +dt*(delta_h*((Lw*Lw)/(2.0*(Lw-WM)))*((psi_1_n0.dx(0)**2)*M11+dot(hat_psi_star.dx(0), (2.0*MN1*psi_1_n0.dx(0)\
                                                                                                       +dot(MNN,hat_psi_star.dx(0)))))\
                     -((1.0/(Lw-WM))*(Lw*Lw)*delta_h.dx(0))*( psi_1_n0.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_star)) \
                                                             +dot(hat_psi_star.dx(0), (DN1*psi_1_n0 + dot(DNN, hat_psi_star))))\
                     +(1.0/h_n1)*(delta_h.dx(0)*((1.0/(Lw-WM))*h_n1.dx(0)*(Lw*Lw))\
                                  -(delta_h/h_n1)*( (Lw*Lw)*(h_n1.dx(0)**2)/(2.0*(Lw-WM))))*(psi_1_n0*psi_1_n0*S11 \
                                                                                             + 2.0*dot(hat_psi_star,SN1)*psi_1_n0\
                                                                                             +dot(hat_psi_star,dot(SNN,hat_psi_star)))\
                     -(0.5*delta_h*(Lw-WM)*H0*H0/(h_n1**2))*(psi_1_n0*psi_1_n0*A11 + 2.0*dot(hat_psi_star,AN1)*psi_1_n0 \
                                                             + dot(hat_psi_star,dot(ANN,hat_psi_star)))\
                     +H0*g*(Lw-WM)*delta_h*(h_n1-H) - H0*psi_1_n0*(x_coord-Lw)*dWM_dt*delta_h.dx(0)))*dx - dt*(Lw*dWM_dt*delta_h*(psi_1_n0*I1 + dot(hat_psi_star,IN)))*ds(1) - dt*(delta_h*(psi_1_n0*G1+dot(hat_psi_star,GN))*hu_fe*Lw/h_n1)*ds(2)

    VP = LinearVariationalProblem(A_psi_s, L_psi_s, psi_1_n1)
    return VP



def energy(t, beach_dw, E_file, Lw, H0, g, WM, dWM, h_n0, psi_1_n0, hat_psi_n0, A11, A1N, AN1, ANN, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, S1N, SN1, SNN, I1, IN, I3):
    
    E = assemble(( 0.5*(h_n0/(Lw-WM))*(Lw**2)*((psi_1_n0.dx(0)**2)*M11 \
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
                  +0.5*H0*g*(Lw-WM)*(h_n0-beach_dw)**2)*dx)
                  
    E_file.write('%-10s %-10s\n' % (str(t), str(E)))
