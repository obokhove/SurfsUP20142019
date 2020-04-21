

from firedrake import *

"""
    ************************************************************************************************************************
    *                                  Weak Formulations for the symplectic Euler scheme                                   *
    ************************************************************************************************************************ """
#--------------------------------------------------------------------------------------------------------------------------#
#                         Step 1 : Update h at time t^{n+1} and psi_i at time t^* simulataneously:                         #
#__________________________________________________________________________________________________________________________#

def WF_h_SE(dim, n_z, g, H, H0, Lw, WM, dWM_dy, dWM_dt, dt, delta_psi, delta_hat_star, h_n0, h_n1, x_coord, psi_1_n0, hat_psi_star, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN):
    # Cf. definitions in JCP article
    Ww = Lw - WM
    Xx = x_coord - Lw
    Uu = Xx*dWM_dy
    Vv = Lw*Lw + Uu**2
    if dim == "3D":
        WF_h = (H0*delta_psi*(h_n1-h_n0)*Ww/dt \
                -(h_n1*(Vv/Ww)*(psi_1_n0.dx(0)*M11 + dot(hat_psi_star.dx(0),MN1))*delta_psi.dx(0)\
                  +Ww*h_n1*(psi_1_n0.dx(1)*M11+dot(hat_psi_star.dx(1),MN1))*delta_psi.dx(1) \
                  +Uu*h_n1*(delta_psi.dx(0)*(M11*psi_1_n0.dx(1) + dot(M1N,hat_psi_star.dx(1))) \
                                             + delta_psi.dx(1)*(M11*psi_1_n0.dx(0) + dot(MN1, hat_psi_star.dx(0))))\
                  -( (Vv/Ww)*h_n1.dx(0) + Uu*h_n1.dx(1) )*( delta_psi.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_star)) \
                                                +delta_psi*(psi_1_n0.dx(0)*D11 + dot(hat_psi_star.dx(0),DN1)))\
                  -( Ww*h_n1.dx(1) + Uu*h_n1.dx(0))*(delta_psi.dx(1)*(D11*psi_1_n0 + dot(D1N,hat_psi_star))\
                                                +delta_psi*(psi_1_n0.dx(1)*D11 + dot(hat_psi_star.dx(1),DN1)))\
                  +(1/h_n1)*( (Vv/Ww)*h_n1.dx(0)**2 + Ww*(h_n1.dx(1)**2  \
                                                + 2.0*Uu*h_n1.dx(0)*h_n1.dx(1) )*(psi_1_n0*S11 + dot(hat_psi_star,SN1))*delta_psi\
                  +(Ww*H0*H0/h_n1)*(psi_1_n0*A11 + dot(hat_psi_star,AN1))*delta_psi -delta_psi*H0*Xx*dWM_dt*h_n1.dx(0)))*dx \
                    - (delta_psi*Lw*dWM_dt*h_n1*I1)*ds(1)
                    
        WF_hat_psi_star= (h_n1*(Vv/Ww)*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_n0.dx(0)+dot(MNN,hat_psi_star.dx(0))))\
                          +Ww*h_n1*elem_mult(delta_hat_star.dx(1),(MN1*psi_1_n0.dx(1)+dot(MNN,hat_psi_star.dx(1))))\
                          +Uu*h_n1*(elem_mult((dot(hat_psi_star.dx(0),MNN)+psi_1_n0.dx(0)*M1N),delta_hat_star.dx(1))\
                                    +elem_mult(delta_hat_star.dx(0),(MN1*psi_1_n0.dx(1)+dot(MNN,hat_psi_star.dx(1)))))\
                          -( (Vv/Ww)*h_n1.dx(0) + Uu*h_n1.dx(1) )*(elem_mult(delta_hat_star,(psi_1_n0.dx(0)*D1N+dot(DNN.T,hat_psi_star.dx(0)))) \
                                    +elem_mult(delta_hat_star.dx(0),(DN1*psi_1_n0+dot(DNN,hat_psi_star)))) \
                          -(Ww*h_n1.dx(1) + Uu*h_n1.dx(0))*(elem_mult(delta_hat_star,(D1N*psi_1_n0.dx(1) + dot(DNN.T,hat_psi_star.dx(1))))\
                                    +elem_mult(delta_hat_star.dx(1),(DN1*psi_1_n0 + dot(DNN,hat_psi_star))))\
                          +(1.0/h_n1)*( (Vv/Ww)*(h_n1.dx(0)**2) + Ww*(h_n1.dx(1)**2) \
                                        +2.0*Uu*h_n1.dx(0)*h_n1.dx(1) )*elem_mult(delta_hat_star,(SN1*psi_1_n0+ dot(SNN,hat_psi_star)))\
                          + (Ww*H0*H0/h_n1)*elem_mult(delta_hat_star,(AN1*psi_1_n0+dot(ANN,hat_psi_star))))

    elif dim=="2D":
        WF_h = (H0*delta_psi*(h_n1-h_n0)*Ww/dt \
                -((h_n1/Ww)*(Lw*Lw)*(psi_1_n0.dx(0)*M11 + dot(hat_psi_star.dx(0),MN1))*delta_psi.dx(0)\
                  -( (1/Ww)*(Lw*Lw)*h_n1.dx(0))*( delta_psi.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_star)) \
                                                      +delta_psi*(psi_1_n0.dx(0)*D11 + dot(hat_psi_star.dx(0),DN1)))\
                  +(1/h_n1)*((1/Ww)*(Lw*Lw)*(h_n1.dx(0)**2))*(psi_1_n0*S11 + dot(hat_psi_star,SN1))*delta_psi\
                  +(Ww*H0*H0/h_n1)*(psi_1_n0*A11 + dot(hat_psi_star,AN1))*delta_psi \
                  -delta_psi*H0*Xx*dWM_dt*h_n1.dx(0)))*dx - (delta_psi*Lw*dWM_dt*h_n1*I1)*ds(1)
                    
        WF_hat_psi_star= ((h_n1/Ww)*(Lw*Lw)*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_n0.dx(0)+dot(MNN,hat_psi_star.dx(0))))\
                          -((Lw*Lw)*h_n1.dx(0)/Ww)*(elem_mult(delta_hat_star, (psi_1_n0.dx(0)*D1N+ dot(DNN.T,hat_psi_star.dx(0)))) \
                                                         + elem_mult(delta_hat_star.dx(0),(DN1*psi_1_n0+dot(DNN,hat_psi_star)))) \
                          +(1.0/h_n1)*((Lw*Lw)*(h_n1.dx(0)**2)/Ww)*elem_mult(delta_hat_star,(SN1*psi_1_n0+ dot(SNN,hat_psi_star)))\
                          + (Ww*H0*H0/h_n1)*elem_mult(delta_hat_star,(AN1*psi_1_n0+dot(ANN,hat_psi_star))))
    
              
    WF_hat_BC = (Lw*dWM_dt*h_n1*elem_mult(delta_hat_star,IN))
    WF_h_psi = WF_h + sum((WF_hat_psi_star[ind])*dx for ind in range(0,n_z)) + sum((WF_hat_BC[ind])*ds(1) for ind in range(0,n_z))

    return WF_h_psi



#----------------------------------------------------------------------------------------------------------------------#
#                                        Step 2 : Update psi_1 at time t^{n+1}:                                        #
#______________________________________________________________________________________________________________________#

def WF_psi_SE(dim, g, H, H0, Lw, WM, WM_n1, dWM_dy, dWM_dt, dt, x_coord, delta_h, psi_1, psi_1_n0, hat_psi_star, h_n1, M11, MN1, MNN, D11, D1N, DN1, DNN,S11, SN1, SNN, A11, AN1, ANN, I1, IN):
    # Cf. definitions in JCP article
    Ww = Lw - WM
    Xx = x_coord - Lw
    Uu = Xx*dWM_dy
    Vv = Lw*Lw + Uu**2
    if dim=="3D":
        A_psi_s = (H0*delta_h*(Lw-WM_n1)*psi_1)*dx
        
        L_psi_s = -(-H0*delta_h*Ww*psi_1_n0 \
                    +dt*( delta_h*( 0.5*(Vv/Ww)*((psi_1_n0.dx(0)**2)*M11 + dot(hat_psi_star.dx(0), (2.0*MN1*psi_1_n0.dx(0)\
                                +dot(MNN,hat_psi_star.dx(0)))))\
                                +0.5*Ww*( (psi_1_n0.dx(1)**2)*M11 + dot(hat_psi_star.dx(1), \
                                                                               (2.0*MN1*psi_1_n0.dx(1) + dot(MNN,hat_psi_star.dx(1)))))\
                                +Uu*( psi_1_n0.dx(0)*(M11*psi_1_n0.dx(1) + dot(MN1,hat_psi_star.dx(1))) \
                                + dot(hat_psi_star.dx(0), (MN1*psi_1_n0.dx(1) + dot(MNN,hat_psi_star.dx(1))))))\
                         -((Vv/Ww)*delta_h.dx(0) + Uu*delta_h.dx(1))*( psi_1_n0.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_star)) \
                                +dot(hat_psi_star.dx(0), (DN1*psi_1_n0 + dot(DNN, hat_psi_star))))\
                         -(Ww*delta_h.dx(1) + Uu*delta_h.dx(0))*( psi_1_n0.dx(1)*(D11*psi_1_n0 + dot(D1N,hat_psi_star))\
                                +dot(hat_psi_star.dx(1), (DN1*psi_1_n0 + dot(DNN, hat_psi_star))))\
                         +(1.0/h_n1)*(delta_h.dx(0)*( (Vv/Ww)*h_n1.dx(0) + Uu*h_n1.dx(1) )\
                                -(delta_h/h_n1)*( 0.5*(Vv/Ww)*(h_n1.dx(0)**2) + 0.5*Ww*(h_n1.dx(1)**2) + Uu*h_n1.dx(0)*h_n1.dx(1) )\
                                + delta_h.dx(1)*( Ww*h_n1.dx(1) + Uu*h_n1.dx(0) ))*( psi_1_n0*psi_1_n0*S11 + 2.0*dot(hat_psi_star,SN1)*psi_1_n0\
                                                                                     + dot(hat_psi_star,dot(SNN,hat_psi_star)) )\
                         -(0.5*delta_h*Ww*H0*H0/(h_n1**2))*(psi_1_n0*psi_1_n0*A11 + 2.0*dot(hat_psi_star,AN1)*psi_1_n0 \
                                + dot(hat_psi_star,dot(ANN,hat_psi_star)))\
                         +H0*g*Ww*delta_h*(h_n1-H) - H0*psi_1_n0*Xx*dWM_dt*delta_h.dx(0)))*dx - dt*(Lw*dWM_dt*delta_h*(psi_1_n0*I1 + dot(hat_psi_star,IN)))*ds(1)

    elif dim=="2D":
        
        A_psi_s = (H0*delta_h*(Lw-WM_n1)*psi_1)*dx
        
        L_psi_s = -(-H0*delta_h*Ww*psi_1_n0 \
                    +dt*(delta_h*((Lw*Lw)/(2.0*Ww))*((psi_1_n0.dx(0)**2)*M11+dot(hat_psi_star.dx(0), (2.0*MN1*psi_1_n0.dx(0)\
                                                                                                           +dot(MNN,hat_psi_star.dx(0)))))\
                         -((1.0/Ww)*(Lw*Lw)*delta_h.dx(0))*( psi_1_n0.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_star)) \
                                                                 +dot(hat_psi_star.dx(0), (DN1*psi_1_n0 + dot(DNN, hat_psi_star))))\
                         +(1.0/h_n1)*(delta_h.dx(0)*((1.0/Ww)*h_n1.dx(0)*(Lw*Lw))\
                                    - (delta_h/h_n1)*( (Lw*Lw)*(h_n1.dx(0)**2)/(2.0*Ww)))*(psi_1_n0*psi_1_n0*S11 \
                                    + 2.0*dot(hat_psi_star,SN1)*psi_1_n0 + dot(hat_psi_star,dot(SNN,hat_psi_star)))\
                         -(0.5*delta_h*Ww*H0*H0/(h_n1**2))*(psi_1_n0*psi_1_n0*A11 + 2.0*dot(hat_psi_star,AN1)*psi_1_n0 \
                                    + dot(hat_psi_star,dot(ANN,hat_psi_star)))\
                         +H0*g*Ww*delta_h*(h_n1-H) - H0*psi_1_n0*Xx*dWM_dt*delta_h.dx(0)))*dx - dt*(Lw*dWM_dt*delta_h*(psi_1_n0*I1 + dot(hat_psi_star,IN)))*ds(1)

    return A_psi_s, L_psi_s


#----------------------------------------------------------------------------------------------------------------------#
#                                        Step 3 : Update psi_i at time t^{n+1}:                                        #
#______________________________________________________________________________________________________________________#
def WF_hat_psi_SE(dim, H, H0, g, n_z, Lw, x_coord, WM, dWM_dt, dWM_dy, dt, delta_hat_psi, hat_psi, h_n0, psi_1_n0, M11, MN1, MNN, D11, D1N, DN1, DNN,S11, SN1, SNN, A11, AN1, ANN, I1, IN):
    # Cf. definitions in JCP article
    Ww = Lw - WM
    Xx = x_coord - Lw
    Uu = Xx*dWM_dy
    Vv = Lw*Lw + Uu**2
    if dim=="3D":
        a_hat_psi =(h_n0*(Vv/Ww)*elem_mult(delta_hat_psi.dx(0), dot(MNN,hat_psi.dx(0)))\
                    +Ww*h_n0*elem_mult(delta_hat_psi.dx(1),dot(MNN,hat_psi.dx(1)))\
                    +Uu*h_n0*( elem_mult(dot(hat_psi.dx(0),MNN),delta_hat_psi.dx(1)) + elem_mult(delta_hat_psi.dx(0),dot(MNN,hat_psi.dx(1))) )\
                    -( (Vv/Ww)*h_n0.dx(0) + Uu*h_n0.dx(1))*( elem_mult(delta_hat_psi,dot(DNN.T,hat_psi.dx(0))) \
                                                             + elem_mult(delta_hat_psi.dx(0),dot(DNN,hat_psi)) ) \
                    -(Ww*h_n0.dx(1) + Uu*h_n0.dx(0))*( elem_mult(delta_hat_psi, dot(DNN.T,hat_psi.dx(1)))\
                                                       + elem_mult(delta_hat_psi.dx(1),dot(DNN,hat_psi)) )\
                    +(1.0/h_n0)*( (Vv/Ww)*(h_n0.dx(0)**2)+Ww*(h_n0.dx(1)**2)+2.0*Uu*h_n0.dx(0)*h_n0.dx(1) )*elem_mult(delta_hat_psi,dot(SNN,hat_psi))\
                    +(Ww*H0*H0/h_n0)*elem_mult(delta_hat_psi,dot(ANN,hat_psi)))
        
        L_hat_psi =-(h_n0*(Vv/Ww)*elem_mult(delta_hat_psi.dx(0), MN1*psi_1_n0.dx(0))\
                     +Ww*h_n0*elem_mult(delta_hat_psi.dx(1), MN1*psi_1_n0.dx(1))\
                     +Uu*h_n0*(elem_mult(delta_hat_psi.dx(1), MN1*psi_1_n0.dx(0))\
                                                +elem_mult(delta_hat_psi.dx(0),MN1*psi_1_n0.dx(1)))\
                     -( (Vv/Ww)*h_n0.dx(0) + Uu*h_n0.dx(1))*( elem_mult(delta_hat_psi, D1N*psi_1_n0.dx(0)) \
                                                              + elem_mult(delta_hat_psi.dx(0),DN1*psi_1_n0) ) \
                     -(Ww*h_n0.dx(1) + Uu*h_n0.dx(0))*( elem_mult(delta_hat_psi, D1N*psi_1_n0.dx(1))\
                                                        +elem_mult(delta_hat_psi.dx(1),DN1*psi_1_n0) )\
                     +(1.0/h_n0)*( (Vv/Ww)*(h_n0.dx(0)**2)+Ww*(h_n0.dx(1)**2)+2.0*Uu*h_n0.dx(0)*h_n0.dx(1) )*elem_mult(delta_hat_psi,SN1*psi_1_n0)\
                     +(Ww*H0*H0/h_n0)*elem_mult(delta_hat_psi,AN1*psi_1_n0))
                     
        hat_psi_BC = -(Lw*dWM_dt*h_n0*elem_mult(delta_hat_psi,IN))
            
    elif dim=="2D":
        a_hat_psi =(h_n0*(Lw*Lw/Ww)*elem_mult(delta_hat_psi.dx(0),dot(MNN,hat_psi.dx(0)))\
                    -(Lw*Lw/Ww)*h_n0.dx(0)*( elem_mult(delta_hat_psi, dot(DNN.T,hat_psi.dx(0)))+elem_mult(delta_hat_psi.dx(0),dot(DNN,hat_psi)) ) \
                    +(1.0/h_n0)*((Lw*Lw)*(h_n0.dx(0)**2)/Ww)*elem_mult(delta_hat_psi,dot(SNN,hat_psi))\
                    +(Ww*H0*H0/h_n0)*elem_mult(delta_hat_psi,dot(ANN,hat_psi)))     
            
        L_hat_psi =-((h_n0/Ww)*(Lw*Lw)*elem_mult(delta_hat_psi.dx(0), MN1*psi_1_n0.dx(0))\
                     -((Lw*Lw)*h_n0.dx(0)/Ww)*(elem_mult(delta_hat_psi, D1N*psi_1_n0.dx(0)) \
                                                    + elem_mult(delta_hat_psi.dx(0),DN1*psi_1_n0)) \
                     +(1.0/h_n0)*( (Lw*Lw)*(h_n0.dx(0)**2)/Ww)*elem_mult(delta_hat_psi,SN1*psi_1_n0)\
                     + (Ww*H0*H0/h_n0)*elem_mult(delta_hat_psi,AN1*psi_1_n0))
        
        hat_psi_BC = -(Lw*dWM_dt*h_n0*elem_mult(delta_hat_psi,IN))


    A_hat = sum((a_hat_psi[ind])*dx for ind in range(0,n_z))
    L_hat = sum((L_hat_psi[ind])*dx for ind in range(0,n_z)) + sum((hat_psi_BC[ind])*ds(1) for ind in range(0,n_z))
    return A_hat, L_hat

"""
    ************************************************************************************************************************
    *                                   Weak Formulations for the Stormer-Verlet scheme                                    *
    ************************************************************************************************************************ """
#----------------------------------------------------------------------------------------------------------------------#
#                                        Step 1 : Update psi_1^{n+1/2} and psi_i^*:                                        #
#______________________________________________________________________________________________________________________#

def WF_psi_half_SV(dim, n_z, g, H, H0, Lw, x_coord, WM, WM_half, dWM_dy, dWM_dt, dWM_half_dy, dWM_half_dt, dt, delta_psi, delta_hat_star, psi_1_n0, psi_1_half, hat_psi_star, h_n0, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN):
    # Cf. definitions in JCP article
    Wwhalf = Lw-WM_half
    Ww = Lw-WM
    Xx = x_coord-Lw
    Uu = Xx*dWM_half_dy
    Vv = Lw*Lw+Uu**2
    if dim=="3D":
        WF_psi_s = ((H0*delta_psi*Wwhalf*(psi_1_half) -delta_psi*Ww*psi_1_n0)/(0.5*dt) \
                    +(delta_psi*( 0.5*(Vv/Wwhalf)*M11*((psi_1_half.dx(0)**2)+dot(hat_psi_star.dx(0),(2.0*MN1*psi_1_half.dx(0)+dot(MNN,hat_psi_star.dx(0)))))\
                                  +0.5*Wwhalf*( (psi_1_half.dx(1)**2)*M11 + dot(hat_psi_star.dx(1),(2.0*MN1*psi_1_half.dx(1) + dot(MNN,hat_psi_star.dx(1)))))\
                                  +Uu*( psi_1_half.dx(0)*(M11*psi_1_half.dx(1) + dot(MN1,hat_psi_star.dx(1))) \
                                        + dot(hat_psi_star.dx(0), (MN1*psi_1_half.dx(1) + dot(MNN,hat_psi_star.dx(1))))))\
                      -((Vv/Wwhalf)**delta_psi.dx(0) + Uu*delta_psi.dx(1))*( psi_1_half.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_star)) \
                                                                             +dot(hat_psi_star.dx(0), (DN1*psi_1_half + dot(DNN, hat_psi_star))))\
                      -(Wwhalf*delta_psi.dx(1) + Uu*delta_psi.dx(0))*( psi_1_half.dx(1)*(D11*psi_1_half + dot(D1N,hat_psi_star))+dot(hat_psi_star.dx(1), (DN1*psi_1_half + dot(DNN, hat_psi_star))))\
                      +(1.0/h_n0)*(delta_psi.dx(0)*((Vv/Wwhalf)*h_n0.dx(0) + Uu*h_n0.dx(1))\
                                   -(delta_psi/h_n0)*( 0.5*(Vv/Wwhalf)*h_n0.dx(0)**2 + 0.5*Wwhalf*(h_n0.dx(1)**2)+Uu*h_n0.dx(0)*h_n0.dx(1) )\
                                   + delta_psi.dx(1)*( Wwhalf*h_n0.dx(1) + Uu*h_n0.dx(0))*(psi_1_half*psi_1_half*S11 + 2.0*dot(hat_psi_star,SN1)*psi_1_half+dot(hat_psi_star,dot(SNN,hat_psi_star)))\
                                   -(0.5*delta_psi*Wwhalf*H0*H0/(h_n0**2))*( psi_1_half*psi_1_half*A11 + 2.0*dot(hat_psi_star,AN1)*psi_1_half+dot(hat_psi_star,dot(ANN,hat_psi_star)) )\
                                   +H0*g*Wwhalf*delta_psi*(h_n0-H) - H0*psi_1_half*Xx*dWM_half_dt*delta_psi.dx(0)))*dx \
                    + (Lw*dWM_half_dt*delta_psi*(psi_1_half*I1 + dot(hat_psi_star,IN)))*ds(1)

        WF_hat_psi_star= (h_n0*(Vv/Wwhalf)*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_half.dx(0)+dot(MNN,hat_psi_star.dx(0))))\
                  + Wwhalf*h_n0*elem_mult(delta_hat_star.dx(1),(MN1*psi_1_half.dx(1)+dot(MNN,hat_psi_star.dx(1))))\
                  + Uu*h_n0*( elem_mult((dot(hat_psi_star.dx(0),MNN)+psi_1_half.dx(0)*M1N),delta_hat_star.dx(1))\
                              + elem_mult(delta_hat_star.dx(0),(MN1*psi_1_half.dx(1)+dot(MNN,hat_psi_star.dx(1)))) )\
                  -( (Vv/Wwhalf)*h_n0.dx(0) + Uu*h_n0.dx(1))*( elem_mult(delta_hat_star,(psi_1_half.dx(0)*D1N+dot(DNN.T,hat_psi_star.dx(0)))) \
                                                               + elem_mult(delta_hat_star.dx(0),(DN1*psi_1_half+dot(DNN,hat_psi_star))) ) \
                  -( Wwhalf*h_n0.dx(1) + Uu*h_n0.dx(0))*( elem_mult(delta_hat_star,(D1N*psi_1_half.dx(1)+ dot(DNN.T,hat_psi_star.dx(1))))\
                                                          + elem_mult(delta_hat_star.dx(1),(DN1*psi_1_half + dot(DNN,hat_psi_star))) )\
                  +(1.0/h_n0)*( (Vv/Wwhalf)*(h_n0.dx(0)**2) + Wwhalf*(h_n0.dx(1)**2) \
                                +2.0*h_n0.dx(0)*h_n0.dx(1)*Uu )*elem_mult(delta_hat_star,(SN1*psi_1_half + dot(SNN,hat_psi_star)))\
                  + (Wwhalf*H0*H0/h_n0)*elem_mult(delta_hat_star,(AN1*psi_1_half+dot(ANN,hat_psi_star))))
    
    if dim=="2D":
        WF_psi_s = ((H0*delta_psi*Wwhalf*(psi_1_half) -delta_psi*Ww*psi_1_n0)/(0.5*dt) \
                    +(delta_psi*( ((Lw*Lw)/(2.0*Wwhalf))*((psi_1_half.dx(0)**2)*M11 +dot(hat_psi_star.dx(0), \
                                                            (2.0*MN1*psi_1_half.dx(0)+dot(MNN,hat_psi_star.dx(0))))))\
                      -((1.0/Wwhalf)*(Lw*Lw)*delta_psi.dx(0))*( psi_1_half.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_star)) \
                                                                +dot(hat_psi_star.dx(0), (DN1*psi_1_half + dot(DNN, hat_psi_star))))\
                      +(1.0/h_n0)*(delta_psi.dx(0)*((1.0/Wwhalf)*h_n0.dx(0)*(Lw*Lw))\
                                   -(delta_psi/h_n0)*( (Lw*Lw)*(h_n0.dx(0)**2)/(2.0*Wwhalf)))*(psi_1_half*psi_1_half*S11 \
                                                        +2.0*dot(hat_psi_star,SN1)*psi_1_half+dot(hat_psi_star,dot(SNN,hat_psi_star)))\
                      -(0.5*delta_psi*Wwhalf*H0*H0/(h_n0**2))*(psi_1_half*psi_1_half*A11 + 2.0*dot(hat_psi_star,AN1)*psi_1_half \
                                                                     + dot(hat_psi_star,dot(ANN,hat_psi_star)))\
                      +H0*g*Wwhalf*delta_psi*(h_n0-H) - H0*psi_1_half*Xx*dWM_half_dt*delta_psi.dx(0)))*dx \
                        + (Lw*dWM_half_dt*delta_psi*(psi_1_half*I1 + dot(hat_psi_star,IN)))*ds(1)                      
                      
        WF_hat_psi_star= ((h_n0/Wwhalf)*(Lw*Lw)*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_half.dx(0)+dot(MNN,hat_psi_star.dx(0))))\
                        -((Lw*Lw)*h_n0.dx(0)/Wwhalf)*( elem_mult(delta_hat_star,(psi_1_half.dx(0)*D1N+dot(DNN.T,hat_psi_star.dx(0)))) \
                                                       + elem_mult(delta_hat_star.dx(0),(DN1*psi_1_half+dot(DNN,hat_psi_star))) ) \
                        +(1.0/h_n0)*((Lw*Lw)*(h_n0.dx(0)**2)/Wwhalf)*elem_mult(delta_hat_star,(SN1*psi_1_half+ dot(SNN,hat_psi_star)))\
                        + (Wwhalf*H0*H0/h_n0)*elem_mult(delta_hat_star,(AN1*psi_1_half+dot(ANN,hat_psi_star))))

    WF_hat_BC_star = (Lw*dWM_half_dt*h_n0*elem_mult(delta_hat_star,IN))
    
    WF_psi_star = WF_psi_s + sum((WF_hat_psi_star[ind])*dx for ind in range(0,n_z)) + sum((WF_hat_BC_star[ind])*ds(1) for ind in range(0,n_z))

    return WF_psi_star


#--------------------------------------------------------------------------------------------------------------------------#
#                       Step 2 : Update h^{n+1} and psi_i at time t^** simulataneously:                       #
#__________________________________________________________________________________________________________________________#

def WF_h_SV(dim, n_z, Lw, H0, g, dt, x_coord, WM, WM_half, dWM_half_dy, dWM_half_dt, delta_psi, delta_hat_star, h_n0, h_n1, psi_1_half, hat_psi_star, hat_psi_aux, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN):
    # Cf. definitions in JCP article
    Wwhalf = Lw-WM_half
    Ww = Lw-WM
    Xx = x_coord-Lw
    Uu = Xx*dWM_half_dy
    Vv = Lw*Lw+Uu**2
    if dim == "3D":
        WF_h = (H0*delta_psi*(h_n1-h_n0)*Wwhalf/dt \
                -0.5*((Vv/Wwhalf)*(h_n0*(psi_1_half.dx(0)*M11 + dot(hat_psi_star.dx(0),MN1)) + h_n1*(psi_1_half.dx(0)*M11 \
                                                                                       + dot(hat_psi_aux.dx(0),MN1)))*delta_psi.dx(0)\
                      +Wwhalf*(h_n0*(psi_1_half.dx(1)*M11+dot(hat_psi_star.dx(1),MN1)) \
                                     +h_n1*(psi_1_half.dx(1)*M11+dot(hat_psi_aux.dx(1),MN1)))*delta_psi.dx(1)
                      +Uu*(h_n0*(delta_psi.dx(0)*(M11*psi_1_half.dx(1) + dot(M1N,hat_psi_star.dx(1))) \
                                                       + delta_psi.dx(1)*(M11*psi_1_half.dx(0) + dot(MN1, hat_psi_star.dx(0))))\
                                                 +h_n1*(delta_psi.dx(0)*(M11*psi_1_half.dx(1) + dot(M1N,hat_psi_aux.dx(1))) \
                                                        + delta_psi.dx(1)*(M11*psi_1_half.dx(0) + dot(MN1, hat_psi_aux.dx(0)))))\
                      -( (Vv/Wwhalf)*h_n0.dx(0) + Uu*h_n0.dx(1))*( delta_psi.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_star)) \
                                                               +delta_psi*(psi_1_half.dx(0)*D11 + dot(hat_psi_star.dx(0),DN1)))\
                      -( (Vv/Wwhalf)*h_n1.dx(0) + Uu*h_n1.dx(1))*( delta_psi.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_aux)) \
                                                               +delta_psi*(psi_1_half.dx(0)*D11 + dot(hat_psi_aux.dx(0),DN1)))\
                      -( Wwhalf*h_n0.dx(1) + Uu*h_n0.dx(0))*(delta_psi.dx(1)*(D11*psi_1_half + dot(D1N,hat_psi_star))\
                                        + delta_psi*(psi_1_half.dx(1)*D11 + dot(hat_psi_star.dx(1),DN1)))\
                      -( Wwhalf*h_n1.dx(1) + Uu*h_n1.dx(0))*(delta_psi.dx(1)*(D11*psi_1_half + dot(D1N,hat_psi_aux))\
                                        + delta_psi*(psi_1_half.dx(1)*D11 + dot(hat_psi_aux.dx(1),DN1)))\
                      +(1/h_n0)*( (Vv/Wwhalf)*(h_n0.dx(0)**2) + Wwhalf*(h_n0.dx(1)**2) \
                                  + 2.0*Uu*h_n0.dx(0)*h_n0.dx(1))*(psi_1_half*S11 + dot(hat_psi_star,SN1))*delta_psi\
                      +(1/h_n1)*( (Vv/Wwhalf)*(h_n1.dx(0)**2) +Wwhalf*(h_n1.dx(1)**2) \
                                 + 2.0*Uu*h_n1.dx(0)*h_n1.dx(1))*(psi_1_half*S11 + dot(hat_psi_aux,SN1))*delta_psi\
                      +(Wwhalf*H0*H0/h_n0)*(psi_1_half*A11 + dot(hat_psi_star,AN1))*delta_psi \
                      +(Wwhalf*H0*H0/h_n1)*(psi_1_half*A11 + dot(hat_psi_aux,AN1))*delta_psi \
                      -delta_psi*H0*Xx*dWM_half_dt*(h_n0.dx(0)+h_n1.dx(0))))*dx \
                    -0.5*(delta_psi*Lw*dWM_half_dt*(h_n0*I1 + h_n1*I1))*ds(1)
                    
        WF_hat_psi_aux= (h_n1*(Vv/Wwhalf)*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_half.dx(0)+dot(MNN,hat_psi_aux.dx(0))))\
                         +Wwhalf*h_n1*elem_mult(delta_hat_star.dx(1),(MN1*psi_1_half.dx(1)+dot(MNN,hat_psi_aux.dx(1))))\
                         +Uu*h_n1*(elem_mult((dot(hat_psi_aux.dx(0),MNN)+psi_1_half.dx(0)*M1N),delta_hat_star.dx(1))\
                                                         +elem_mult(delta_hat_star.dx(0),(MN1*psi_1_half.dx(1)+dot(MNN,hat_psi_aux.dx(1)))))\
                         -( (Vv/Wwhalf)*h_n1.dx(0) + Uu*h_n1.dx(1))*( elem_mult(delta_hat_star,(psi_1_half.dx(0)*D1N+dot(DNN.T,hat_psi_aux.dx(0)))) \
                                                                      + elem_mult(delta_hat_star.dx(0),(DN1*psi_1_half+dot(DNN,hat_psi_aux))) ) \
                         -(Wwhalf*h_n1.dx(1) + Uu*h_n1.dx(0))*(elem_mult(delta_hat_star, (D1N*psi_1_half.dx(1) + dot(DNN.T,hat_psi_aux.dx(1))))\
                                                            + elem_mult(delta_hat_star.dx(1),(DN1*psi_1_half + dot(DNN,hat_psi_aux))))\
                         +(1.0/h_n1)*( (Vv/Wwhalf)*(h_n1.dx(0)**2) + Wwhalf*(h_n1.dx(1)**2) + 2.0*Uu*h_n1.dx(0)*h_n1.dx(1) )*elem_mult(delta_hat_star,(SN1*psi_1_half+dot(SNN,hat_psi_aux)))\
                         + (Wwhalf*H0*H0/h_n1)*elem_mult(delta_hat_star,(AN1*psi_1_half+dot(ANN,hat_psi_aux))))

    elif dim=="2D":
        WF_h = (H0*delta_psi*(h_n1-h_n0)*Wwhalf/dt \
                -0.5*((h_n0/Wwhalf)*(Lw*Lw)*(psi_1_half.dx(0)*M11 + dot(hat_psi_star.dx(0),MN1))*delta_psi.dx(0)\
                      -( (1/Wwhalf)*(Lw*Lw)*h_n0.dx(0))*( delta_psi.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_star)) \
                                                               +delta_psi*(psi_1_half.dx(0)*D11 + dot(hat_psi_star.dx(0),DN1)))\
                      +(1/h_n0)*( (1/Wwhalf)*(Lw*Lw)*(h_n0.dx(0)**2))*(psi_1_half*S11 + dot(hat_psi_star,SN1))*delta_psi\
                      +(Wwhalf*H0*H0/h_n0)*(psi_1_half*A11 + dot(hat_psi_star,AN1))*delta_psi - delta_psi*H0*Xx*dWM_half_dt*h_n0.dx(0))
                -0.5*((h_n1/(Wwhalf)*(Lw*Lw)*(psi_1_half.dx(0)*M11 + dot(hat_psi_aux.dx(0),MN1))*delta_psi.dx(0)\
                      -( (1/Wwhalf)*(Lw*Lw)*h_n1.dx(0))*( delta_psi.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_aux)) \
                                                               +delta_psi*(psi_1_half.dx(0)*D11 + dot(hat_psi_aux.dx(0),DN1)))\
                      +(1/h_n1)*( (1/Wwhalf)*(Lw*Lw)*(h_n1.dx(0)**2))*(psi_1_half*S11 + dot(hat_psi_aux,SN1))*delta_psi\
                       +(Wwhalf*H0*H0/h_n1)*(psi_1_half*A11 + dot(hat_psi_aux,AN1))*delta_psi \
                       -delta_psi*H0*Xx*dWM_half_dt*h_n1.dx(0)) )*dx \
                    -0.5*(delta_psi*Lw*dWM_half_dt*h_n0*I1 + delta_psi*Lw*dWM_half_dt*h_n1*I1)*ds(1)
                        
        WF_hat_psi_aux= ((h_n1/Wwhalf)*(Lw*Lw)*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_half.dx(0)+dot(MNN,hat_psi_aux.dx(0))))\
                         -((Lw*Lw)*h_n1.dx(0)/Wwhalf)*(elem_mult(delta_hat_star, (psi_1_half.dx(0)*D1N+dot(DNN.T,hat_psi_aux.dx(0)))) \
                                                             + elem_mult(delta_hat_star.dx(0),(DN1*psi_1_half+dot(DNN,hat_psi_aux)))) \
                         +(1.0/h_n1)*( (Lw*Lw)*(h_n1.dx(0)**2)/Wwhalf)*elem_mult(delta_hat_star,(SN1*psi_1_half+ dot(SNN,hat_psi_aux)))\
                         + (Wwhalf*H0*H0/h_n1)*elem_mult(delta_hat_star,(AN1*psi_1_half+dot(ANN,hat_psi_aux))))


    WF_hat_BC_aux = (Lw*dWM_half_dt*h_n1*elem_mult(delta_hat_star,IN))
    WF_h_psi = WF_h + sum((WF_hat_psi_aux[ind])*dx for ind in range(0,n_z)) + sum((WF_hat_BC_aux[ind])*ds(1) for ind in range(0,n_z))

    return WF_h_psi

#----------------------------------------------------------------------------------------------------------------------#
#                                        Step 3 : Update psi_1 at time t^{n+1}:                                        #
#______________________________________________________________________________________________________________________#

def WF_psi_n1_SV(dim, H0, H, g, x_coord, delta_h, Lw, WM_n1, WM_half, dt, psi_1_half, psi_1, dWM_half_dt, dWM_half_dy, hat_psi_aux, h_n1, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN ):
    # Cf. definitions in JCP article
    Wwhalf = Lw - WM_half
    Xx = x_coord - Lw
    Uu = Xx*dWM_half_dy
    Vv = Lw*Lw + Uu**2
    if dim=="3D":
        a_psi_1 = (H0*delta_h*(Lw-WM_n1)*(psi_1)/(0.5*dt))*dx
        L_psi_1 = -( -H0*delta_h*Wwhalf*psi_1_half/(0.5*dt) \
                    +(delta_h*(0.5*(Vv/Wwhalf)*((psi_1_half.dx(0)**2)*M11 + dot(hat_psi_aux.dx(0), (2.0*MN1*psi_1_half.dx(0) + dot(MNN,hat_psi_aux.dx(0)))))\
                               +0.5*Wwhalf*( (psi_1_half.dx(1)**2)*M11 + dot(hat_psi_aux.dx(1), (2.0*MN1*psi_1_half.dx(1) + dot(MNN,hat_psi_aux.dx(1)))))\
                               +Uu*( psi_1_half.dx(0)*(M11*psi_1_half.dx(1) + dot(MN1,hat_psi_aux.dx(1))) \
                                     + dot(hat_psi_aux.dx(0), (MN1*psi_1_half.dx(1) + dot(MNN,hat_psi_aux.dx(1))))))\
                      -((Vv/Wwhalf)*delta_h.dx(0) + Uu*delta_h.dx(1))*( psi_1_half.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_aux)) \
                                                                        + dot(hat_psi_aux.dx(0), (DN1*psi_1_half + dot(DNN, hat_psi_aux))) )\
                      -(Wwhalf*delta_h.dx(1) + Uu*delta_h.dx(0))*( psi_1_half.dx(1)*(D11*psi_1_half + dot(D1N,hat_psi_aux))+dot(hat_psi_aux.dx(1), (DN1*psi_1_half + dot(DNN, hat_psi_aux))) )\
                      +(1.0/h_n1)*(delta_h.dx(0)*((Vv/Wwhalf)*h_n1.dx(0) + Uu*h_n1.dx(1))\
                                   -(delta_h/h_n1)*( 0.5*(Vv/Wwhalf)*(h_n1.dx(0)**2) + 0.5*Wwhalf*(h_n1.dx(1)**2) + Uu*h_n1.dx(0)*h_n1.dx(1) )\
                                   + delta_h.dx(1)*( Wwhalf*h_n1.dx(1) + Uu*h_n1.dx(0))*( psi_1_half*psi_1_half*S11 \
                                                                                          + 2.0*dot(hat_psi_aux,SN1)*psi_1_half + dot(hat_psi_aux,dot(SNN,hat_psi_aux)) )\
                      -(0.5*delta_h*Wwhalf*H0*H0/(h_n1**2))*( psi_1_half*psi_1_half*A11 + 2.0*dot(hat_psi_aux,AN1)*psi_1_half \
                                                              + dot(hat_psi_aux,dot(ANN,hat_psi_aux)) )\
                                   +H0*g*Wwhalf*delta_h*(h_n1-H) - H0*psi_1_half*Xx*dWM_half_dt*delta_h.dx(0)) )*dx \
                        - (Lw*dWM_half_dt*delta_h*(psi_1_half*I1 + dot(hat_psi_aux,IN)))*ds(1)

    elif dim=="2D":
        a_psi_1 = (H0*delta_h*(Lw-WM_n1)*(psi_1)/(0.5*dt))*dx
        L_psi_1 = -( -H0*delta_h*Wwhalf*psi_1_half/(0.5*dt) \
                    +(delta_h*(((Lw*Lw)/(2.0*Wwhalf))*( (psi_1_half.dx(0)**2)*M11 +dot(hat_psi_aux.dx(0),(2.0*MN1*psi_1_half.dx(0) + dot(MNN,hat_psi_aux.dx(0))))) )\
                      -((1.0/Wwhalf)*(Lw*Lw)*delta_h.dx(0))*( psi_1_half.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_aux)) \
                                                              +dot(hat_psi_aux.dx(0), (DN1*psi_1_half + dot(DNN, hat_psi_aux))) )\
                      +(1.0/h_n1)*(delta_h.dx(0)*((1.0/Wwhalf)*h_n1.dx(0)*(Lw*Lw))\
                                   -(delta_h/h_n1)*( (Lw*Lw)*(h_n1.dx(0)**2)/(2.0*Wwhalf)))*(psi_1_half*psi_1_half*S11 \
                                                     +2.0*dot(hat_psi_aux,SN1)*psi_1_half+dot(hat_psi_aux,dot(SNN,hat_psi_aux)) )\
                      -(0.5*delta_h*Wwhalf*H0*H0/(h_n1**2))*( psi_1_half*psi_1_half*A11 + 2.0*dot(hat_psi_aux,AN1)*psi_1_half \
                                                              + dot(hat_psi_aux,dot(ANN,hat_psi_aux)))\
                      +H0*g*Wwhalf*delta_h*(h_n1-H) - H0*psi_1_half*(x_coord-Lw)*dWM_half_dt*delta_h.dx(0)) )*dx \
                        - (Lw*dWM_half_dt*delta_h*(psi_1_half*I1 + dot(hat_psi_aux,IN)))*ds(1)

    return a_psi_1, L_psi_1


#----------------------------------------------------------------------------------------------------------------------#
#                                        Step 4 : Update psi_i at time t^{n+1}:                                        #
#______________________________________________________________________________________________________________________#

def WF_hat_psi_SV(dim, n_z, Lw, H0, H, WM, x_coord, dt, dWM_dt, dWM_dy, delta_hat_psi, hat_psi, h_n0, psi_1_n0, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN ):
    Ww = Lw-WM
    Xx = x_coord-Lw
    Uu = Xx*dWM_dy
    Vv = Lw*Lw+Uu**2
    if dim=="3D":
        a_hat_psi =( h_n0*(Vv/Ww)*elem_mult(delta_hat_psi.dx(0), dot(MNN,hat_psi.dx(0)))\
                     +Ww*h_n0*elem_mult(delta_hat_psi.dx(1),dot(MNN,hat_psi.dx(1)))\
                     +Uu*h_n0*(elem_mult(dot(hat_psi.dx(0),MNN),delta_hat_psi.dx(1))\
                                               +elem_mult(delta_hat_psi.dx(0),dot(MNN,hat_psi.dx(1))))\
                     -( (Vv/Ww)*h_n0.dx(0) + Uu*h_n0.dx(1))*( elem_mult(delta_hat_psi, dot(DNN.T,hat_psi.dx(0))) \
                                                              + elem_mult(delta_hat_psi.dx(0),dot(DNN,hat_psi)) ) \
                     -(Ww*h_n0.dx(1) + Uu*h_n0.dx(0))*( elem_mult(delta_hat_psi, dot(DNN.T,hat_psi.dx(1)))\
                                                        +elem_mult(delta_hat_psi.dx(1),dot(DNN,hat_psi)) )\
                     +(1.0/h_n0)*( (Vv/Ww)*(h_n0.dx(0)**2) + Ww*(h_n0.dx(1)**2) + 2.0*Uu*h_n0.dx(0)*h_n0.dx(1) )*elem_mult(delta_hat_psi,dot(SNN,hat_psi))\
                     + (Ww*H0*H0/h_n0)*elem_mult(delta_hat_psi,dot(ANN,hat_psi)))
        
        L_hat_psi =-( h_n0*(Vv/Ww)*elem_mult(delta_hat_psi.dx(0), MN1*psi_1_n0.dx(0)) + Ww*h_n0*elem_mult(delta_hat_psi.dx(1), MN1*psi_1_n0.dx(1))\
                      +Uu*h_n0*( elem_mult(delta_hat_psi.dx(1), MN1*psi_1_n0.dx(0)) + elem_mult(delta_hat_psi.dx(0),MN1*psi_1_n0.dx(1)) )\
                      -( (Vv/Ww)*h_n0.dx(0) + Uu*h_n0.dx(1))*( elem_mult(delta_hat_psi, D1N*psi_1_n0.dx(0)) + elem_mult(delta_hat_psi.dx(0),DN1*psi_1_n0) )\
                      -( Ww*h_n0.dx(1) + Uu*h_n0.dx(0))*(elem_mult(delta_hat_psi, D1N*psi_1_n0.dx(1)) + elem_mult(delta_hat_psi.dx(1),DN1*psi_1_n0))\
                      +(1.0/h_n0)*( (Vv/Ww)*h_n0.dx(0)**2 + Ww*h_n0.dx(1)**2 + 2.0*Uu*h_n0.dx(0)*h_n0.dx(1) )*elem_mult(delta_hat_psi,SN1*psi_1_n0)\
                      +(Ww*H0**2/h_n0)*elem_mult(delta_hat_psi,AN1*psi_1_n0))
            
    elif dim=="2D":
         a_hat_psi =((h_n0/Ww)*Lw**2*elem_mult(delta_hat_psi.dx(0), dot(MNN,hat_psi.dx(0)))\
                     -(Lw**2*h_n0.dx(0)/Ww)*( elem_mult(delta_hat_psi, dot(DNN.T,hat_psi.dx(0))) + elem_mult(delta_hat_psi.dx(0),dot(DNN,hat_psi)) )\
                     +(1.0/h_n0)*((Lw*Lw)*(h_n0.dx(0)**2)/Ww)*elem_mult(delta_hat_psi,dot(SNN,hat_psi))\
                     + (Ww*H0**2/h_n0)*elem_mult(delta_hat_psi,dot(ANN,hat_psi)))

         L_hat_psi =-((h_n0/Ww)*Lw**2*elem_mult(delta_hat_psi.dx(0), MN1*psi_1_n0.dx(0))\
                      -(Lw**2*h_n0.dx(0)/Ww)*( elem_mult(delta_hat_psi, D1N*psi_1_n0.dx(0)) + elem_mult(delta_hat_psi.dx(0),DN1*psi_1_n0) ) \
                      +(1.0/h_n0)*( Lw**2*(h_n0.dx(0)**2)/Ww)*elem_mult(delta_hat_psi,SN1*psi_1_n0) + (Ww*H0**2/h_n0)*elem_mult(delta_hat_psi,AN1*psi_1_n0) )

    hat_psi_BC = -(Lw*dWM_dt*h_n0*elem_mult(delta_hat_psi,IN))
    A_hat = sum((a_hat_psi[ind])*dx for ind in range(0,n_z))
    L_hat = sum((L_hat_psi[ind])*dx for ind in range(0,n_z)) + sum((hat_psi_BC[ind])*ds(1) for ind in range(0,n_z))

    return A_hat, L_hat


