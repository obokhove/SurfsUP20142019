

from firedrake import *

def op(V,W,U,f,g):
    op_grad = (V/W)*f.dx(0)*g.dx(0)+W*f.dx(1)*g.dx(1)+U*(f.dx(0)*g.dx(1)+f.dx(1)*g.dx(0))
    return op_grad

def op_vec(V,W,U,f,g,A):
#    op_grad = (V/W)*dot(A,g.dx(0))*f.dx(0)+W*dot(A,g.dx(1))*f.dx(1)+U*(dot(A,g.dx(0))*f.dx(1)+dot(A,g.dx(1))*f.dx(0))
    op_grad = (V/W)*dot(g.dx(0),A)*f.dx(0)+W*dot(g.dx(1),A)*f.dx(1)+U*(dot(g.dx(0),A)*f.dx(1)+dot(g.dx(1),A)*f.dx(0))
    return op_grad

#def op_elem(V,W,U,f,g,A):
#    f = delta_hat_star
#    g = hat_psi_star
#    op_grad =(V/W)*elem_mult(f.dx(0),dot(MNN,g.dx(0)))+W*elem_mult(f.dx(1),dot(MNN,g.dx(1)))+U*(elem_mult(dot(f.dx(1),MNN),g.dx(0))\
#                 +elem_mult(f.dx(0),dot(MNN,g.dx(1))))\
#    -( (V/W)*h_n1.dx(0)+U*h_n1.dx(1))*(elem_mult(delta_hat_star,dot(DNN.T,hat_psi_star.dx(0)) ) \
#                                   + elem_mult(delta_hat_star.dx(0),dot(DNN,hat_psi_star))) \
#    -(W*h_n1.dx(1) + U*h_n1.dx(0))*(elem_mult(delta_hat_star,(dot(DNN.T,hat_psi_star.dx(1))))\
#                                    +elem_mult(delta_hat_star.dx(1),(dot(DNN,hat_psi_star))))\

"""
    ************************************************************************************************************************
    *                                  Weak Formulations for the symplectic Euler scheme                                   *
    ************************************************************************************************************************ """
#--------------------------------------------------------------------------------------------------------------------------#
#                         Step 1 : Update h at time t^{n+1} and psi_i at time t^* simulataneously:                         #
#__________________________________________________________________________________________________________________________#

def WF_h_SE(dim, n_z, g, b, H0, Lw, WM, dWM_dy, dWM_dt, dt, delta_psi, delta_hat_star, h_n0, h_n1, x_coord, psi_1_n0, hat_psi_star, A11, AN1, ANN, B11, B1N, BN1, BNN, C11, CN1, CNN, D11, D1N, DN1, DNN, M11,M1N, MN1, MNN, S11, SN1, SNN, I1, IN):
    X = x_coord-Lw
    U = X*dWM_dy
    W=(Lw-WM)
    V=Lw*Lw+U**2
    
    if dim == "3D":
        WF_h = (H0*delta_psi*(h_n1-h_n0)*W/dt \
                -(h_n1*( M11*op(V,W,U,delta_psi,psi_1_n0)+ op_vec(V,W,U,delta_psi,hat_psi_star,MN1))\
                  -op(V,W,U,h_n1,delta_psi)*(D11*psi_1_n0 + dot(D1N,hat_psi_star))\
                  -H0*op(V,W,U,b,delta_psi)*(B11*psi_1_n0 + dot(B1N,hat_psi_star))\
                  - delta_psi*(D11*op(V,W,U,psi_1_n0,h_n1) + H0*B11*op(V,W,U,psi_1_n0,b)\
                               +op_vec(V,W,U,h_n1,hat_psi_star,DN1)+H0*op_vec(V,W,U,b,hat_psi_star,BN1))
                  +(1/h_n1)*op(V,W,U,h_n1,h_n1)*(psi_1_n0*S11 + dot(hat_psi_star,SN1))*delta_psi \
                  +(H0*H0/h_n1)*(op(V,W,U,b,b)+W)*(psi_1_n0*A11 + dot(hat_psi_star,AN1))*delta_psi \
                  +(2*H0/h_n1)*op(V,W,U,b,h_n1)*(psi_1_n0*C11 + dot(hat_psi_star,CN1))*delta_psi \
                  -delta_psi*H0*X*dWM_dt*h_n1.dx(0)))*dx \
                    - (delta_psi*Lw*dWM_dt*h_n1*I1)*ds(1)

#        WF_hat_psi_star= (h_n1*(elem_mult(op(delta_hat_star, psi_1_n0),MN1) + op_elem(V,W,U, delta_hat_star, hat_psi_star,MNN)) \
#                          -op(V,W,U,h_n1,psi_1_n0)*elem_mult(delta_hat_star,D1N) - H0*op(V,W,U,b,psi_1_n0)*elem_mult(delta_hat_star,B1N)\
#                          - (elem_mult(op(V,W,U,h_n1,delta_hat_star),DN1) +H0*elem_mult(op(V,W,U,b,delta_hat_star),BN1))*psi_1_n0\
#                          -elem_mult(delta_hat_star,op_vec(V,W,U,h_n1,hat_psi_star,DNN.T))\
#                          -H0*elem_mult(delta_hat_star,op_vec(V,W,U,b,hat_psi_star,BNN.T))\
#                          -elem_mult(op(V,W,U,h_n1,delta_hat_star), dot(DNN,hat_psi_star))\
#                          -H0*elem_mult(op(V,W,U,b,delta_hat_star), dot(BNN,hat_psi_star))\
#                          +(1.0/h_n1)*op(V,W,U,h_n1,h_n1)*elem_mult(delta_hat_star,(SN1*psi_1_n0+ dot(SNN,hat_psi_star)))\
#                          + (H0*H0/h_n1)*(W+op(V,W,U,b,b))*elem_mult(delta_hat_star,(AN1*psi_1_n0+dot(ANN,hat_psi_star)))\
#                          + (2*H0/h_n1)*(W+op(V,W,U,b,h_n1))*elem_mult(delta_hat_star,(CN1*psi_1_n0+dot(CNN,hat_psi_star))))
        WF_hat_psi_star= (h_n1*(V/W)*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_n0.dx(0)+dot(MNN,hat_psi_star.dx(0))))\
                          +W*h_n1*elem_mult(delta_hat_star.dx(1),(MN1*psi_1_n0.dx(1)+dot(MNN,hat_psi_star.dx(1))))\
                          +U*h_n1*(elem_mult((dot(hat_psi_star.dx(0),MNN)+psi_1_n0.dx(0)*M1N),delta_hat_star.dx(1))\
                                                     +elem_mult(delta_hat_star.dx(0),(MN1*psi_1_n0.dx(1)+dot(MNN,hat_psi_star.dx(1)))))\
                          -((V/W)*h_n1.dx(0)+U*h_n1.dx(1))*(elem_mult(delta_hat_star, \
                                                                       (psi_1_n0.dx(0)*D1N+dot(DNN.T,hat_psi_star.dx(0)))) \
                                                              + elem_mult(delta_hat_star.dx(0),(DN1*psi_1_n0+dot(DNN,hat_psi_star)))) \
                          -(W*h_n1.dx(1) + U*h_n1.dx(0))*(elem_mult(delta_hat_star,(D1N*psi_1_n0.dx(1)+dot(DNN.T,hat_psi_star.dx(1))))\
                                                          +elem_mult(delta_hat_star.dx(1),(DN1*psi_1_n0 + dot(DNN,hat_psi_star))))\
                          -H0*((V/W)*b.dx(0)+U*b.dx(1))*(elem_mult(delta_hat_star,(psi_1_n0.dx(0)*B1N+dot(BNN.T,hat_psi_star.dx(0)))) \
                                                            + elem_mult(delta_hat_star.dx(0),(BN1*psi_1_n0+dot(BNN,hat_psi_star)))) \
                          -H0*(W*b.dx(1) + U*b.dx(0))*(elem_mult(delta_hat_star,(B1N*psi_1_n0.dx(1)+dot(BNN.T,hat_psi_star.dx(1))))\
                                                          +elem_mult(delta_hat_star.dx(1),(BN1*psi_1_n0 + dot(BNN,hat_psi_star))))\
                          +(1.0/h_n1)*((V/W)*(h_n1.dx(0)**2) +W*(h_n1.dx(1)**2) +2.0*h_n1.dx(0)*h_n1.dx(1)*U)*elem_mult(delta_hat_star,(SN1*psi_1_n0\
                                                                                                                 + dot(SNN,hat_psi_star)))\
                          + (H0*H0/h_n1)*(W+op(V,W,U,b,b))*elem_mult(delta_hat_star,(AN1*psi_1_n0+dot(ANN,hat_psi_star)))\
                          + (2*H0/h_n1)*(op(V,W,U,b,h_n1))*elem_mult(delta_hat_star,(CN1*psi_1_n0+dot(CNN,hat_psi_star))))

    elif dim=="2D":
        WF_h = (H0*delta_psi*(h_n1-h_n0)*W/dt \
                -((h_n1/W)*(Lw*Lw)*(psi_1_n0.dx(0)*M11 + dot(hat_psi_star.dx(0),MN1))*delta_psi.dx(0)\
                  -( (1/W)*(Lw*Lw)*h_n1.dx(0))*( delta_psi.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_star)) \
                                                      +delta_psi*(psi_1_n0.dx(0)*D11 + dot(hat_psi_star.dx(0),DN1)))\
                  +(1/h_n1)*((1/W)*(Lw*Lw)*(h_n1.dx(0)**2))*(psi_1_n0*S11 + dot(hat_psi_star,SN1))*delta_psi\
                  +(W*H0*H0/h_n1)*(psi_1_n0*A11 + dot(hat_psi_star,AN1))*delta_psi \
                  -delta_psi*H0*X*dWM_dt*h_n1.dx(0)))*dx - (delta_psi*Lw*dWM_dt*h_n1*I1)*ds(1)
                    
        WF_hat_psi_star= ((h_n1/W)*(Lw*Lw)*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_n0.dx(0)+dot(MNN,hat_psi_star.dx(0))))\
                          -((Lw*Lw)*h_n1.dx(0)/W)*(elem_mult(delta_hat_star, (psi_1_n0.dx(0)*D1N+ dot(DNN.T,hat_psi_star.dx(0)))) \
                                                         + elem_mult(delta_hat_star.dx(0),(DN1*psi_1_n0+dot(DNN,hat_psi_star)))) \
                          +(1.0/h_n1)*((Lw*Lw)*(h_n1.dx(0)**2)/W)*elem_mult(delta_hat_star,(SN1*psi_1_n0+ dot(SNN,hat_psi_star)))\
                          + (W*H0*H0/h_n1)*elem_mult(delta_hat_star,(AN1*psi_1_n0+dot(ANN,hat_psi_star))))
    
              
    WF_hat_BC = (Lw*dWM_dt*h_n1*elem_mult(delta_hat_star,IN))
    WF_h_psi = WF_h + sum((WF_hat_psi_star[ind])*dx for ind in range(0,n_z)) + sum((WF_hat_BC[ind])*ds(1) for ind in range(0,n_z))

    return WF_h_psi



#----------------------------------------------------------------------------------------------------------------------#
#                                        Step 2 : Update psi_1 at time t^{n+1}:                                        #
#______________________________________________________________________________________________________________________#

def WF_psi_SE(dim, g, H, b, H0, Lw, WM, WM_n1, dWM_dy, dWM_dt, dt, x_coord, delta_h, psi_1, psi_1_n0, hat_psi_star, h_n1, A11, AN1, ANN, B11, B1N, BN1, BNN, C11, CN1, CNN, D11, D1N, DN1, DNN, M11, MN1, MNN, S11, SN1, SNN, I1, IN):
    X = x_coord-Lw
    U = X*dWM_dy
    W=(Lw-WM)
    V=Lw*Lw+U**2
    if dim=="3D":
#        A_psi_s = (H0*delta_h*(Lw-WM_n1)*psi_1)*dx
#        
#        L_psi_s = -(-H0*delta_h*(Lw-WM)*psi_1_n0 \
#                    +dt*(delta_h*( ((Lw*Lw+((x_coord-Lw)*dWM_dy)**2)/(2.0*(Lw-WM)))*((psi_1_n0.dx(0)**2)*M11 \
#                                                                                     +dot(hat_psi_star.dx(0), (2.0*MN1*psi_1_n0.dx(0)\
#                                                                                                               +dot(MNN,hat_psi_star.dx(0)))))\
#                                  +0.5*(Lw-WM)*( (psi_1_n0.dx(1)**2)*M11 + dot(hat_psi_star.dx(1), \
#                                                                               (2.0*MN1*psi_1_n0.dx(1) + dot(MNN,hat_psi_star.dx(1)))))\
#                                  +(x_coord-Lw)*dWM_dy*( psi_1_n0.dx(0)*(M11*psi_1_n0.dx(1) + dot(MN1,hat_psi_star.dx(1))) \
#                                                        + dot(hat_psi_star.dx(0), (MN1*psi_1_n0.dx(1) + dot(MNN,hat_psi_star.dx(1))))))\
#                         -((1.0/(Lw-WM))*(Lw*Lw+((x_coord-Lw)*dWM_dy)**2)*delta_h.dx(0) \
#                           + (x_coord-Lw)*dWM_dy*delta_h.dx(1))*( psi_1_n0.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_star)) \
#                                                                 +dot(hat_psi_star.dx(0), (DN1*psi_1_n0 + dot(DNN, hat_psi_star))))\
#                         -((Lw-WM)*delta_h.dx(1) + (x_coord-Lw)*delta_h.dx(0)*dWM_dy)*( psi_1_n0.dx(1)*(D11*psi_1_n0 + dot(D1N,hat_psi_star))\
#                                                                                       +dot(hat_psi_star.dx(1), (DN1*psi_1_n0 + dot(DNN, hat_psi_star))))\
#                         +(1.0/h_n1)*(delta_h.dx(0)*((1.0/(Lw-WM))*h_n1.dx(0)*(Lw*Lw+((x_coord-Lw)*dWM_dy)**2) + h_n1.dx(1)*(x_coord-Lw)*dWM_dy)\
#                                      -(delta_h/h_n1)*( (Lw*Lw+((x_coord-Lw)*dWM_dy)**2)*(h_n1.dx(0)**2)/(2.0*(Lw-WM)) + 0.5*(Lw-WM)*(h_n1.dx(1)**2)\
#                                                       + h_n1.dx(0)*h_n1.dx(1)*(x_coord-Lw)*dWM_dy )\
#                                      + delta_h.dx(1)*( (Lw-WM)*h_n1.dx(1) + h_n1.dx(0)*(x_coord-Lw)*dWM_dy))*(psi_1_n0*psi_1_n0*S11 \
#                                                                                                               + 2.0*dot(hat_psi_star,SN1)*psi_1_n0\
#                                                                                                               +dot(hat_psi_star,dot(SNN,hat_psi_star)))\
#                         -(0.5*delta_h*(Lw-WM)*H0*H0/(h_n1**2))*(psi_1_n0*psi_1_n0*A11 + 2.0*dot(hat_psi_star,AN1)*psi_1_n0 \
#                                                                 + dot(hat_psi_star,dot(ANN,hat_psi_star)))\
#                         +H0*g*(Lw-WM)*delta_h*(h_n1-H) - H0*psi_1_n0*(x_coord-Lw)*dWM_dt*delta_h.dx(0)))*dx - dt*(Lw*dWM_dt*delta_h*(psi_1_n0*I1 + dot(hat_psi_star,IN)))*ds(1)
        A_psi_s = (H0*delta_h*(Lw-WM_n1)*psi_1)*dx

        L_psi_s = -(-H0*delta_h*W*psi_1_n0 \
                    +dt*( 0.5*delta_h*(M11*op(V,W,U,psi_1_n0,psi_1_n0) + 2*op_vec(V,W,U,psi_1_n0,hat_psi_star,MN1)\
                                       +(V/W)*(dot(hat_psi_star.dx(0), (dot(MNN,hat_psi_star.dx(0)))))\
                                       +W*dot(hat_psi_star.dx(1), (dot(MNN,hat_psi_star.dx(1))))\
                                       +2*U*psi_1_n0.dx(0)*dot(hat_psi_star.dx(0), dot(MNN,hat_psi_star.dx(1))))\
                         -op(V,W,U,psi_1_n0,delta_h)*(D11*psi_1_n0 + dot(D1N,hat_psi_star))\
                         -op_vec(V,W,U,delta_h,hat_psi_star,DN1)*psi_1_n0 -dot(op_vec(V,W,U,delta_h,hat_psi_star,DNN),hat_psi_star)\
                         +(1.0/h_n1)*(op(V,W,U,delta_h,h_n1)\
                                      -(delta_h/h_n1)*0.5*op(V,W,U,h_n1,h_n1))*(psi_1_n0*psi_1_n0*S11 \
                                                                                + 2.0*dot(hat_psi_star,SN1)*psi_1_n0\
                                                                                +dot(hat_psi_star,dot(SNN,hat_psi_star)))\
                         +(H0/h_n1)*(op(V,W,U,delta_h,b)\
                                      -(delta_h/h_n1)*op(V,W,U,h_n1,b))*(psi_1_n0*psi_1_n0*C11 + 2.0*dot(hat_psi_star,CN1)*psi_1_n0\
                                                                                +dot(hat_psi_star,dot(CNN,hat_psi_star)))\
                         -(0.5*delta_h*H0*H0/(h_n1**2))*(W+op(V,W,U,b,b))*(psi_1_n0*psi_1_n0*A11 + 2.0*dot(hat_psi_star,AN1)*psi_1_n0 \
                                                                 + dot(hat_psi_star,dot(ANN,hat_psi_star)))\
                         +H0*g*W*delta_h*(h_n1-H) - H0*psi_1_n0*X*dWM_dt*delta_h.dx(0)))*dx - dt*(Lw*dWM_dt*delta_h*(psi_1_n0*I1 + dot(hat_psi_star,IN)))*ds(1)

    elif dim=="2D":
        
        A_psi_s = (H0*delta_h*(Lw-WM_n1)*psi_1)*dx
        
        L_psi_s = -(-H0*delta_h*W*psi_1_n0 \
                    +dt*(delta_h*((Lw*Lw)/(2.0*W))*((psi_1_n0.dx(0)**2)*M11+dot(hat_psi_star.dx(0), (2.0*MN1*psi_1_n0.dx(0)\
                                                                                                           +dot(MNN,hat_psi_star.dx(0)))))\
                         -((1.0/W)*(Lw*Lw)*delta_h.dx(0))*( psi_1_n0.dx(0)*(D11*psi_1_n0 + dot(D1N,hat_psi_star)) \
                                                                 +dot(hat_psi_star.dx(0), (DN1*psi_1_n0 + dot(DNN, hat_psi_star))))\
                         +(1.0/h_n1)*(delta_h.dx(0)*((1.0/W)*h_n1.dx(0)*(Lw*Lw))\
                                      -(delta_h/h_n1)*( (Lw*Lw)*(h_n1.dx(0)**2)/(2.0*W)))*(psi_1_n0*psi_1_n0*S11 \
                                                                                                 + 2.0*dot(hat_psi_star,SN1)*psi_1_n0\
                                                                                                 +dot(hat_psi_star,dot(SNN,hat_psi_star)))\
                         -(0.5*delta_h*W*H0*H0/(h_n1**2))*(psi_1_n0*psi_1_n0*A11 + 2.0*dot(hat_psi_star,AN1)*psi_1_n0 \
                                                                 + dot(hat_psi_star,dot(ANN,hat_psi_star)))\
                         +H0*g*W*delta_h*(h_n1-H) - H0*psi_1_n0*X*dWM_dt*delta_h.dx(0)))*dx - dt*(Lw*dWM_dt*delta_h*(psi_1_n0*I1 + dot(hat_psi_star,IN)))*ds(1)


    return A_psi_s, L_psi_s


#----------------------------------------------------------------------------------------------------------------------#
#                                        Step 3 : Update psi_i at time t^{n+1}:                                        #
#______________________________________________________________________________________________________________________#
def WF_hat_psi_SE(dim, b, H, H0, g, n_z, Lw, x_coord, WM, dWM_dt, dWM_dy, dt, delta_hat_psi, hat_psi, h_n0, psi_1_n0, A11, AN1, ANN, B11, B1N, BN1, BNN, C11, C1N, CN1, CNN, D11, D1N, DN1, DNN, M11, MN1, MNN, S11, SN1, SNN, I1, IN):
    X = x_coord-Lw
    U = X*dWM_dy
    W=(Lw-WM)
    V=Lw*Lw+U**2
    if dim=="3D":
#        a_hat_psi =((h_n0/(Lw-WM))*(Lw*Lw+((x_coord-Lw)*dWM_dy)**2)*elem_mult(delta_hat_psi.dx(0), dot(MNN,hat_psi.dx(0)))\
#                    +(Lw-WM)*h_n0*elem_mult(delta_hat_psi.dx(1),dot(MNN,hat_psi.dx(1)))\
#                    +(x_coord-Lw)*dWM_dy*h_n0*(elem_mult(dot(hat_psi.dx(0),MNN),delta_hat_psi.dx(1))\
#                                               +elem_mult(delta_hat_psi.dx(0),dot(MNN,hat_psi.dx(1))))\
#                    -( (Lw*Lw+((x_coord-Lw)*(dWM_dy))**2)*h_n0.dx(0)/(Lw-WM) \
#                      +(x_coord-Lw)*dWM_dy*h_n0.dx(1))*(elem_mult(delta_hat_psi, dot(DNN.T,hat_psi.dx(0))) \
#                                                        + elem_mult(delta_hat_psi.dx(0),dot(DNN,hat_psi))) \
#                    -((Lw-WM)*h_n0.dx(1) + (x_coord-Lw)*dWM_dy*h_n0.dx(0))*(elem_mult(delta_hat_psi, dot(DNN.T,hat_psi.dx(1)))\
#                                                                            +elem_mult(delta_hat_psi.dx(1),dot(DNN,hat_psi)))\
#                    +(1.0/h_n0)*( (Lw*Lw+((x_coord-Lw)*(dWM_dy))**2)*(h_n0.dx(0)**2)/(Lw-WM) +(Lw-WM)*(h_n0.dx(1)**2) \
#                                 +2.0*h_n0.dx(0)*h_n0.dx(1)*(x_coord-Lw)*dWM_dy)*elem_mult(delta_hat_psi,dot(SNN,hat_psi))\
#                    + ((Lw-WM)*H0*H0/h_n0)*elem_mult(delta_hat_psi,dot(ANN,hat_psi)))

        a_hat_psi =(h_n0*((V/W)*elem_mult(delta_hat_psi.dx(0), dot(MNN,hat_psi.dx(0)))\
                    +W*elem_mult(delta_hat_psi.dx(1),dot(MNN,hat_psi.dx(1)))\
                    +U*(elem_mult(dot(hat_psi.dx(0),MNN),delta_hat_psi.dx(1))\
                                               +elem_mult(delta_hat_psi.dx(0),dot(MNN,hat_psi.dx(1)))))\
                    -((V/W)*h_n0.dx(0)+U*h_n0.dx(1))*(elem_mult(delta_hat_psi, dot(DNN.T,hat_psi.dx(0))) \
                                                        + elem_mult(delta_hat_psi.dx(0),dot(DNN,hat_psi))) \
                    -(W*h_n0.dx(1) + U*h_n0.dx(0))*(elem_mult(delta_hat_psi, dot(DNN.T,hat_psi.dx(1)))\
                                                                            +elem_mult(delta_hat_psi.dx(1),dot(DNN,hat_psi)))\
                    -H0*((V/W)*b.dx(0)+U*b.dx(1))*(elem_mult(delta_hat_psi, dot(BNN.T,hat_psi.dx(0))) \
                                                        + elem_mult(delta_hat_psi.dx(0),dot(BNN,hat_psi))) \
                    -H0*(W*b.dx(1) + U*b.dx(0))*(elem_mult(delta_hat_psi, dot(BNN.T,hat_psi.dx(1)))\
                                                                            +elem_mult(delta_hat_psi.dx(1),dot(BNN,hat_psi)))\
                    +(1.0/h_n0)*((V/W)*(h_n0.dx(0)**2) +W*(h_n0.dx(1)**2)+2.0*h_n0.dx(0)*h_n0.dx(1)*U)*elem_mult(delta_hat_psi,dot(SNN,hat_psi))\
                    + (H0*H0/h_n0)*(W+op(V,W,U,b,b))*elem_mult(delta_hat_psi,dot(ANN,hat_psi))\
                    +(2*H0/h_n0)*op(V,W,U,b,h_n0)*elem_mult(delta_hat_psi,dot(CNN,hat_psi)))

        
        L_hat_psi =-((h_n0/W)*V*elem_mult(delta_hat_psi.dx(0), MN1*psi_1_n0.dx(0))\
                     +W*h_n0*elem_mult(delta_hat_psi.dx(1), MN1*psi_1_n0.dx(1))\
                     +U*h_n0*(elem_mult(delta_hat_psi.dx(1), MN1*psi_1_n0.dx(0))\
                                                +elem_mult(delta_hat_psi.dx(0),MN1*psi_1_n0.dx(1)))\
                     -( h_n0.dx(0)*(V/W) +U*h_n0.dx(1))*(elem_mult(delta_hat_psi, D1N*psi_1_n0.dx(0)) \
                                                         + elem_mult(delta_hat_psi.dx(0),DN1*psi_1_n0)) \
                     -(W*h_n0.dx(1) + U*h_n0.dx(0))*(elem_mult(delta_hat_psi, D1N*psi_1_n0.dx(1))\
                                                                             +elem_mult(delta_hat_psi.dx(1),DN1*psi_1_n0))\
                     -H0*(b.dx(0)*(V/W)+U*b.dx(1))*(elem_mult(delta_hat_psi, B1N*psi_1_n0.dx(0)) \
                                                         + elem_mult(delta_hat_psi.dx(0),BN1*psi_1_n0)) \
                     -H0*(W*b.dx(1) + U*b.dx(0))*(elem_mult(delta_hat_psi, B1N*psi_1_n0.dx(1))\
                                                                             +elem_mult(delta_hat_psi.dx(1),BN1*psi_1_n0))\
                     +(1.0/h_n0)*( (V/W)*(h_n0.dx(0)**2) +W*(h_n0.dx(1)**2)+2.0*h_n0.dx(0)*h_n0.dx(1)*U)*elem_mult(delta_hat_psi,SN1*psi_1_n0)\
                     + (H0*H0/h_n0)*(W+op(V,W,U,b,b))*elem_mult(delta_hat_psi,AN1*psi_1_n0)\
                     + (2*H0/h_n0)*op(V,W,U,b,h_n0)*elem_mult(delta_hat_psi,CN1*psi_1_n0))
                     
        hat_psi_BC = -(Lw*dWM_dt*h_n0*elem_mult(delta_hat_psi,IN))
            
    elif dim=="2D":
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

def WF_psi_half_SV(b,dim, n_z, g, H, H0, Lw, x_coord, WM, WM_half, dWM_dy, dWM_dt, dWM_half_dy, dWM_half_dt, dt, delta_psi, delta_hat_star, psi_1_n0, psi_1_half, hat_psi_star, h_n0, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN, B11, B1N, BN1, BNN, C11, C1N, CN1, CNN ):

    if dim=="3D":
        X = x_coord-Lw
        U = X*dWM_half_dy
        W=(Lw-WM)
        W_half=(Lw-WM_half)
        V=Lw*Lw+U**2
        

        WF_psi_s = ((H0*delta_psi*W_half*(psi_1_half) -delta_psi*W*psi_1_n0)/(0.5*dt) \
                    +(delta_psi*( (V/(2.0*W_half))*((psi_1_half.dx(0)**2)*M11 +dot(hat_psi_star.dx(0), (2.0*MN1*psi_1_half.dx(0)\
                                                                                                        +dot(MNN,hat_psi_star.dx(0)))))\
                                 +0.5*W_half*( (psi_1_half.dx(1)**2)*M11 + dot(hat_psi_star.dx(1), \
                                                                               (2.0*MN1*psi_1_half.dx(1) + dot(MNN,hat_psi_star.dx(1)))))\
                                 +X*dWM_half_dy*( psi_1_half.dx(0)*(M11*psi_1_half.dx(1) + dot(MN1,hat_psi_star.dx(1))) \
                                                            + dot(hat_psi_star.dx(0), (MN1*psi_1_half.dx(1) + dot(MNN,hat_psi_star.dx(1))))))\
                      -((1.0/W_half)*V*delta_psi.dx(0)+ U*delta_psi.dx(1))*( psi_1_half.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_star)) \
                                                                                        +dot(hat_psi_star.dx(0), (DN1*psi_1_half + dot(DNN, hat_psi_star))))\
                      -(W_half*delta_psi.dx(1) + U*delta_psi.dx(0))*( psi_1_half.dx(1)*(D11*psi_1_half + dot(D1N,hat_psi_star))\
                                                                                 +dot(hat_psi_star.dx(1), (DN1*psi_1_half + dot(DNN, hat_psi_star))))\
                      +(1.0/h_n0)*(delta_psi.dx(0)*((1.0/W_half)*h_n0.dx(0)*V + h_n0.dx(1)*U)\
                                   -(delta_psi/h_n0)*( V*(h_n0.dx(0)**2)/(2.0*W_half) + 0.5*W_half*(h_n0.dx(1)**2)\
                                                      + h_n0.dx(0)*h_n0.dx(1)*U )\
                                   + delta_psi.dx(1)*( W_half*h_n0.dx(1) + h_n0.dx(0)*U))*(psi_1_half*psi_1_half*S11 + 2.0*dot(hat_psi_star,SN1)*psi_1_half+dot(hat_psi_star,dot(SNN,hat_psi_star)))\
                      +(H0/h_n0)*(delta_psi.dx(0)*((1.0/W_half)*b.dx(0)*V + b.dx(1)*U)\
                                   -(delta_psi/h_n0)*( V*(h_n0.dx(0)*b.dx(0))/(W_half) + W_half*(h_n0.dx(1)*b.dx(1))\
                                                      + b.dx(0)*h_n0.dx(1)*U+ h_n0.dx(0)*b.dx(1)*U )\
                                   + delta_psi.dx(1)*( W_half*b.dx(1) + b.dx(0)*U))*(psi_1_half*psi_1_half*C11 + 2.0*dot(hat_psi_star,CN1)*psi_1_half+dot(hat_psi_star,dot(CNN,hat_psi_star)))\
                      -(0.5*delta_psi*H0*H0/(h_n0**2))*(W_half+op(V,W_half,U,b,b))*(psi_1_half*psi_1_half*A11 + 2.0*dot(hat_psi_star,AN1)*psi_1_half \
                                                                     + dot(hat_psi_star,dot(ANN,hat_psi_star)))\
                      +H0*g*W_half*delta_psi*(h_n0-H) - H0*psi_1_half*X*dWM_half_dt*delta_psi.dx(0)))*dx \
                        + (Lw*dWM_half_dt*delta_psi*(psi_1_half*I1 + dot(hat_psi_star,IN)))*ds(1)
                    
                    
                    
        WF_hat_psi_star= ((h_n0/W_half)*V*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_half.dx(0)+dot(MNN,hat_psi_star.dx(0))))\
                          +W_half*h_n0*elem_mult(delta_hat_star.dx(1),(MN1*psi_1_half.dx(1)+dot(MNN,hat_psi_star.dx(1))))\
                          +X*dWM_half_dy*h_n0*(elem_mult((dot(hat_psi_star.dx(0),MNN)+psi_1_half.dx(0)*M1N),delta_hat_star.dx(1))\
                                                          +elem_mult(delta_hat_star.dx(0),(MN1*psi_1_half.dx(1)+dot(MNN,hat_psi_star.dx(1)))))\
                          -( V*h_n0.dx(0)/W_half +U*h_n0.dx(1))*(elem_mult(delta_hat_star, (psi_1_half.dx(0)*D1N+dot(DNN.T,hat_psi_star.dx(0)))) \
                                                                 + elem_mult(delta_hat_star.dx(0),(DN1*psi_1_half+dot(DNN,hat_psi_star)))) \
                          -(W_half*h_n0.dx(1) + U*h_n0.dx(0))*(elem_mult(delta_hat_star, (D1N*psi_1_half.dx(1)+dot(DNN.T,hat_psi_star.dx(1))))\
                                                               +elem_mult(delta_hat_star.dx(1),(DN1*psi_1_half + dot(DNN,hat_psi_star))))\
                          -H0*(V*b.dx(0)/W_half +U*b.dx(1))*(elem_mult(delta_hat_star, (psi_1_half.dx(0)*B1N+dot(BNN.T,hat_psi_star.dx(0)))) \
                                                             + elem_mult(delta_hat_star.dx(0),(BN1*psi_1_half+dot(BNN,hat_psi_star)))) \
                          -H0*(W_half*b.dx(1) + U*b.dx(0))*(elem_mult(delta_hat_star, (B1N*psi_1_half.dx(1)+dot(BNN.T,hat_psi_star.dx(1))))\
                                                            +elem_mult(delta_hat_star.dx(1),(BN1*psi_1_half + dot(BNN,hat_psi_star))))\
                          +(1.0/h_n0)*op(V,W_half,U,h_n0,h_n0)*elem_mult(delta_hat_star,(SN1*psi_1_half+ dot(SNN,hat_psi_star)))\
                          +(2*H0/h_n0)*op(V,W_half,U,h_n0,b)*elem_mult(delta_hat_star,(CN1*psi_1_half+ dot(CNN,hat_psi_star)))\
                          +(H0*H0/h_n0)*(W_half+op(V,W_half,U,b,b))*elem_mult(delta_hat_star,(AN1*psi_1_half+dot(ANN,hat_psi_star))))
    if dim=="2D":
        WF_psi_s = ((H0*delta_psi*(Lw-WM_half)*(psi_1_half) -delta_psi*W*psi_1_n0)/(0.5*dt) \
                    +(delta_psi*( ((Lw*Lw)/(2.0*(Lw-WM_half)))*((psi_1_half.dx(0)**2)*M11 +dot(hat_psi_star.dx(0), \
                                                                                               (2.0*MN1*psi_1_half.dx(0)\
                                                                                                +dot(MNN,hat_psi_star.dx(0))))))\
                      -((1.0/(Lw-WM_half))*(Lw*Lw)*delta_psi.dx(0))*( psi_1_half.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_star)) \
                                                                     +dot(hat_psi_star.dx(0), (DN1*psi_1_half + dot(DNN, hat_psi_star))))\
                      +(1.0/h_n0)*(delta_psi.dx(0)*((1.0/(Lw-WM_half))*h_n0.dx(0)*(Lw*Lw))\
                                   -(delta_psi/h_n0)*( (Lw*Lw)*(h_n0.dx(0)**2)/(2.0*(Lw-WM_half))))*(psi_1_half*psi_1_half*S11 \
                                                                                                     + 2.0*dot(hat_psi_star,SN1)*psi_1_half\
                                                                                                     +dot(hat_psi_star,dot(SNN,hat_psi_star)))\
                      -(0.5*delta_psi*(Lw-WM_half)*H0*H0/(h_n0**2))*(psi_1_half*psi_1_half*A11 + 2.0*dot(hat_psi_star,AN1)*psi_1_half \
                                                                     + dot(hat_psi_star,dot(ANN,hat_psi_star)))\
                      +H0*g*(Lw-WM_half)*delta_psi*(h_n0-H) - H0*psi_1_half*(x_coord-Lw)*dWM_half_dt*delta_psi.dx(0)))*dx \
                        + (Lw*dWM_half_dt*delta_psi*(psi_1_half*I1 + dot(hat_psi_star,IN)))*ds(1)
                      
                      
                      
        WF_hat_psi_star= ((h_n0/(Lw-WM_half))*(Lw*Lw)*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_half.dx(0)+dot(MNN,hat_psi_star.dx(0))))\
                        -((Lw*Lw)*h_n0.dx(0)/(Lw-WM_half))*(elem_mult(delta_hat_star, \
                                                                      (psi_1_half.dx(0)*D1N+dot(DNN.T,hat_psi_star.dx(0)))) \
                                                            + elem_mult(delta_hat_star.dx(0),(DN1*psi_1_half+dot(DNN,hat_psi_star)))) \
                        +(1.0/h_n0)*((Lw*Lw)*(h_n0.dx(0)**2)/(Lw-WM_half))*elem_mult(delta_hat_star,(SN1*psi_1_half+ dot(SNN,hat_psi_star)))\
                        + ((Lw-WM_half)*H0*H0/h_n0)*elem_mult(delta_hat_star,(AN1*psi_1_half+dot(ANN,hat_psi_star))))

    WF_hat_BC_star = (Lw*dWM_half_dt*h_n0*elem_mult(delta_hat_star,IN))
    
    WF_psi_star = WF_psi_s + sum((WF_hat_psi_star[ind])*dx for ind in range(0,n_z)) + sum((WF_hat_BC_star[ind])*ds(1) for ind in range(0,n_z))

    return WF_psi_star


#--------------------------------------------------------------------------------------------------------------------------#
#                       Step 2 : Update h^{n+1} and psi_i at time t^** simulataneously:                       #
#__________________________________________________________________________________________________________________________#

def WF_h_SV(b,dim, n_z, Lw, H0, g, dt, x_coord, WM, WM_half, dWM_half_dy, dWM_half_dt, delta_psi, delta_hat_star, h_n0, h_n1, psi_1_half, hat_psi_star, hat_psi_aux, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN, B11, B1N, BN1, BNN, C11, C1N, CN1, CNN ):

    if dim == "3D":
        X = x_coord-Lw
        U = X*dWM_half_dy
        W=(Lw-WM)
        W_half=(Lw-WM_half)
        V=Lw*Lw+U**2
        
        WF_h = (H0*delta_psi*(h_n1-h_n0)*W_half/dt \
                -0.5*((1.0/W_half)*V*(h_n0*(psi_1_half.dx(0)*M11 + dot(hat_psi_star.dx(0),MN1))+h_n1*(psi_1_half.dx(0)*M11 \
                                                                                                      + dot(hat_psi_aux.dx(0),MN1)))*delta_psi.dx(0)\
                      +W_half*(h_n0*(psi_1_half.dx(1)*M11+dot(hat_psi_star.dx(1),MN1)) \
                               +h_n1*(psi_1_half.dx(1)*M11+dot(hat_psi_aux.dx(1),MN1)))*delta_psi.dx(1)
                      +U*(h_n0*(delta_psi.dx(0)*(M11*psi_1_half.dx(1) + dot(M1N,hat_psi_star.dx(1))) \
                                + delta_psi.dx(1)*(M11*psi_1_half.dx(0) + dot(MN1, hat_psi_star.dx(0))))\
                          +h_n1*(delta_psi.dx(0)*(M11*psi_1_half.dx(1) + dot(M1N,hat_psi_aux.dx(1))) \
                                 + delta_psi.dx(1)*(M11*psi_1_half.dx(0) + dot(MN1, hat_psi_aux.dx(0)))))\
                      -( (1.0/W_half)*V*h_n0.dx(0)+U*h_n0.dx(1))*( delta_psi.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_star)) \
                                                                  +delta_psi*(psi_1_half.dx(0)*D11 + dot(hat_psi_star.dx(0),DN1)))\
                      -( (1.0/W_half)*V*h_n1.dx(0)+U*h_n1.dx(1))*( delta_psi.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_aux)) \
                                                                  +delta_psi*(psi_1_half.dx(0)*D11 + dot(hat_psi_aux.dx(0),DN1)))\
                      -( W_half*h_n0.dx(1) + U*h_n0.dx(0))*(delta_psi.dx(1)*(D11*psi_1_half + dot(D1N,hat_psi_star))\
                                                            +delta_psi*(psi_1_half.dx(1)*D11+ dot(hat_psi_star.dx(1),DN1)))\
                      -( W_half*h_n1.dx(1) + U*h_n1.dx(0))*(delta_psi.dx(1)*(D11*psi_1_half + dot(D1N,hat_psi_aux))\
                                                            +delta_psi*(psi_1_half.dx(1)*D11 + dot(hat_psi_aux.dx(1),DN1)))\
                      -H0*((1.0/W_half)*V*b.dx(0)+U*b.dx(1))*( delta_psi.dx(0)*(B11*psi_1_half + dot(B1N,hat_psi_star)) \
                                                                  +delta_psi*(psi_1_half.dx(0)*B11 + dot(hat_psi_star.dx(0),BN1)))\
                      -H0*((1.0/W_half)*V*b.dx(0)+U*b.dx(1))*( delta_psi.dx(0)*(B11*psi_1_half + dot(B1N,hat_psi_aux)) \
                                                                  +delta_psi*(psi_1_half.dx(0)*B11 + dot(hat_psi_aux.dx(0),BN1)))\
                      -H0*( W_half*b.dx(1) + U*b.dx(0))*(delta_psi.dx(1)*(B11*psi_1_half + dot(B1N,hat_psi_star))\
                                                            +delta_psi*(psi_1_half.dx(1)*B11+ dot(hat_psi_star.dx(1),BN1)))
                      -H0*( W_half*b.dx(1) + U*b.dx(0))*(delta_psi.dx(1)*(B11*psi_1_half + dot(B1N,hat_psi_aux))\
                                                            +delta_psi*(psi_1_half.dx(1)*B11 + dot(hat_psi_aux.dx(1),BN1)))\
                      +(1/h_n0)*op(V,W_half,U,h_n0,h_n0)*(psi_1_half*S11 + dot(hat_psi_star,SN1))*delta_psi\
                      +(1/h_n1)*op(V,W_half,U,h_n1,h_n1)*(psi_1_half*S11 + dot(hat_psi_aux,SN1))*delta_psi\
                      +(2*H0/h_n0)*op(V,W_half,U,b,h_n0)*(psi_1_half*C11 + dot(hat_psi_star,CN1))*delta_psi\
                      +(2*H0/h_n1)*op(V,W_half,U,b,h_n1)*(psi_1_half*C11 + dot(hat_psi_aux,CN1))*delta_psi\
                      +(H0*H0/h_n0)*(W_half+op(V,W_half,U,b,b))*(psi_1_half*A11 + dot(hat_psi_star,AN1))*delta_psi \
                      +(H0*H0/h_n1)*(W_half+op(V,W_half,U,b,b))*(psi_1_half*A11 + dot(hat_psi_aux,AN1))*delta_psi \
                      -delta_psi*H0*X*dWM_half_dt*(h_n0.dx(0)+h_n1.dx(0))))*dx \
                    -0.5*(delta_psi*Lw*dWM_half_dt*(h_n0*I1 + h_n1*I1))*ds(1)
                    
        WF_hat_psi_aux= ((h_n1/W_half)*V*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_half.dx(0)+dot(MNN,hat_psi_aux.dx(0))))\
                         +W_half*h_n1*elem_mult(delta_hat_star.dx(1),(MN1*psi_1_half.dx(1)+dot(MNN,hat_psi_aux.dx(1))))\
                         +U*h_n1*(elem_mult((dot(hat_psi_aux.dx(0),MNN)+psi_1_half.dx(0)*M1N),delta_hat_star.dx(1))\
                                  +elem_mult(delta_hat_star.dx(0),(MN1*psi_1_half.dx(1)+dot(MNN,hat_psi_aux.dx(1)))))\
                         -( V*h_n1.dx(0)/W_half+U*h_n1.dx(1))*(elem_mult(delta_hat_star, (psi_1_half.dx(0)*D1N+dot(DNN.T,hat_psi_aux.dx(0)))) \
                                                               + elem_mult(delta_hat_star.dx(0),(DN1*psi_1_half+dot(DNN,hat_psi_aux)))) \
                         -(W_half*h_n1.dx(1) + U*h_n1.dx(0))*(elem_mult(delta_hat_star,(D1N*psi_1_half.dx(1)+dot(DNN.T,hat_psi_aux.dx(1))))\
                                                              +elem_mult(delta_hat_star.dx(1),(DN1*psi_1_half + dot(DNN,hat_psi_aux))))\
                         -H0*( V*b.dx(0)/W_half+U*b.dx(1))*(elem_mult(delta_hat_star, (psi_1_half.dx(0)*B1N+dot(BNN.T,hat_psi_aux.dx(0)))) \
                                                               + elem_mult(delta_hat_star.dx(0),(BN1*psi_1_half+dot(BNN,hat_psi_aux)))) \
                         -H0*(W_half*b.dx(1) + U*b.dx(0))*(elem_mult(delta_hat_star,(B1N*psi_1_half.dx(1)+dot(BNN.T,hat_psi_aux.dx(1))))\
                                                              +elem_mult(delta_hat_star.dx(1),(BN1*psi_1_half + dot(BNN,hat_psi_aux))))\
                         +(1.0/h_n1)*op(V,W_half,U,h_n1,h_n1)*elem_mult(delta_hat_star,(SN1*psi_1_half + dot(SNN,hat_psi_aux)))\
                         +(2.0*H0/h_n1)*op(V,W_half,U,b,h_n1)*elem_mult(delta_hat_star,(CN1*psi_1_half + dot(CNN,hat_psi_aux)))\
                         + (H0*H0/h_n1)*(W_half+op(V,W_half,U,b,b))*elem_mult(delta_hat_star,(AN1*psi_1_half+dot(ANN,hat_psi_aux))))

    elif dim=="2D":
        WF_h = (H0*delta_psi*(h_n1-h_n0)*(Lw-WM_half)/dt \
                -0.5*((h_n0/(Lw-WM_half))*(Lw*Lw)*(psi_1_half.dx(0)*M11 + dot(hat_psi_star.dx(0),MN1))*delta_psi.dx(0)\
                      -( (1/(Lw-WM_half))*(Lw*Lw)*h_n0.dx(0))*( delta_psi.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_star)) \
                                                               +delta_psi*(psi_1_half.dx(0)*D11 + dot(hat_psi_star.dx(0),DN1)))\
                      +(1/h_n0)*( (1/(Lw-WM_half))*(Lw*Lw)*(h_n0.dx(0)**2))*(psi_1_half*S11 + dot(hat_psi_star,SN1))*delta_psi\
                      +((Lw-WM_half)*H0*H0/h_n0)*(psi_1_half*A11 + dot(hat_psi_star,AN1))*delta_psi \
                      -delta_psi*H0*(x_coord-Lw)*dWM_half_dt*h_n0.dx(0))
                -0.5*((h_n1/(Lw-WM_half))*(Lw*Lw)*(psi_1_half.dx(0)*M11 + dot(hat_psi_aux.dx(0),MN1))*delta_psi.dx(0)\
                      -( (1/(Lw-WM_half))*(Lw*Lw)*h_n1.dx(0))*( delta_psi.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_aux)) \
                                                               +delta_psi*(psi_1_half.dx(0)*D11 + dot(hat_psi_aux.dx(0),DN1)))\
                      +(1/h_n1)*( (1/(Lw-WM_half))*(Lw*Lw)*(h_n1.dx(0)**2))*(psi_1_half*S11 + dot(hat_psi_aux,SN1))*delta_psi\
                      +((Lw-WM_half)*H0*H0/h_n1)*(psi_1_half*A11 + dot(hat_psi_aux,AN1))*delta_psi \
                      -delta_psi*H0*(x_coord-Lw)*dWM_half_dt*h_n1.dx(0)))*dx \
                    -0.5*(delta_psi*Lw*dWM_half_dt*h_n0*I1 + delta_psi*Lw*dWM_half_dt*h_n1*I1)*ds(1)
                        
        WF_hat_psi_aux= ((h_n1/(Lw-WM_half))*(Lw*Lw)*elem_mult(delta_hat_star.dx(0),(MN1*psi_1_half.dx(0)+dot(MNN,hat_psi_aux.dx(0))))\
                         -((Lw*Lw)*h_n1.dx(0)/(Lw-WM_half))*(elem_mult(delta_hat_star, (psi_1_half.dx(0)*D1N+dot(DNN.T,hat_psi_aux.dx(0)))) \
                                                             + elem_mult(delta_hat_star.dx(0),(DN1*psi_1_half+dot(DNN,hat_psi_aux)))) \
                         +(1.0/h_n1)*( (Lw*Lw)*(h_n1.dx(0)**2)/(Lw-WM_half))*elem_mult(delta_hat_star,(SN1*psi_1_half+ dot(SNN,hat_psi_aux)))\
                         + ((Lw-WM_half)*H0*H0/h_n1)*elem_mult(delta_hat_star,(AN1*psi_1_half+dot(ANN,hat_psi_aux))))


    WF_hat_BC_aux = (Lw*dWM_half_dt*h_n1*elem_mult(delta_hat_star,IN))
    WF_h_psi = WF_h + sum((WF_hat_psi_aux[ind])*dx for ind in range(0,n_z)) + sum((WF_hat_BC_aux[ind])*ds(1) for ind in range(0,n_z))

    return WF_h_psi

#----------------------------------------------------------------------------------------------------------------------#
#                                        Step 3 : Update psi_1 at time t^{n+1}:                                        #
#______________________________________________________________________________________________________________________#

def WF_psi_n1_SV(b,dim, H0, H, g, x_coord, delta_h, Lw, WM_n1, WM_half, dt, psi_1_half, psi_1, dWM_half_dt, dWM_half_dy, hat_psi_aux, h_n1, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN, B11, B1N, BN1, BNN, C11, C1N, CN1, CNN ):
    if dim=="3D":
        X = x_coord-Lw
        U = X*dWM_half_dy
        W_half=(Lw-WM_half)
        V=Lw*Lw+U**2
        
        a_psi_1 = (H0*delta_h*(Lw-WM_n1)*(psi_1)/(0.5*dt))*dx
        L_psi_1 = -( -H0*delta_h*W_half*psi_1_half/(0.5*dt) \
                    +(delta_h*((V/(2.0*W_half))*((psi_1_half.dx(0)**2)*M11 +dot(hat_psi_aux.dx(0), (2.0*MN1*psi_1_half.dx(0)\
                                                                                                    +dot(MNN,hat_psi_aux.dx(0)))))\
                               +0.5*W_half*( (psi_1_half.dx(1)**2)*M11 + dot(hat_psi_aux.dx(1), \
                                                                             (2.0*MN1*psi_1_half.dx(1) + dot(MNN,hat_psi_aux.dx(1)))))\
                               +U*( psi_1_half.dx(0)*(M11*psi_1_half.dx(1) + dot(MN1,hat_psi_aux.dx(1))) \
                                   + dot(hat_psi_aux.dx(0), (MN1*psi_1_half.dx(1) + dot(MNN,hat_psi_aux.dx(1))))))\
                      -((1.0/W_half)*V*delta_h.dx(0) + U*delta_h.dx(1))*( psi_1_half.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_aux)) \
                                                                         +dot(hat_psi_aux.dx(0), (DN1*psi_1_half + dot(DNN, hat_psi_aux))))\
                      -(W_half*delta_h.dx(1) + X*delta_h.dx(0)*dWM_half_dy)*( psi_1_half.dx(1)*(D11*psi_1_half + dot(D1N,hat_psi_aux))\
                                                                             +dot(hat_psi_aux.dx(1), (DN1*psi_1_half + dot(DNN, hat_psi_aux))))\
                      +(1.0/h_n1)*(delta_h.dx(0)*((1.0/W_half)*h_n1.dx(0)*V + h_n1.dx(1)*U)\
                                   -(delta_h/h_n1)*( V*(h_n1.dx(0)**2)/(2.0*W_half) + 0.5*W_half*(h_n1.dx(1)**2)+ h_n1.dx(0)*h_n1.dx(1)*U )\
                                   + delta_h.dx(1)*( W_half*h_n1.dx(1) + h_n1.dx(0)*U))*(psi_1_half*psi_1_half*S11 + 2.0*dot(hat_psi_aux,SN1)*psi_1_half\
                                                                                         +dot(hat_psi_aux,dot(SNN,hat_psi_aux)))\
                      +(H0/h_n1)*(op(V,W_half,U,delta_h,b)-(delta_h/h_n1)*op(V,W_half,U,b,h_n1))*(psi_1_half*psi_1_half*C11 + 2.0*dot(hat_psi_aux,CN1)*psi_1_half\
                                                                                         +dot(hat_psi_aux,dot(CNN,hat_psi_aux)))\
                      -(0.5*delta_h*H0*H0/(h_n1**2))*(W_half+op(V,W_half,U,b,b))*(psi_1_half*psi_1_half*A11 + 2.0*dot(hat_psi_aux,AN1)*psi_1_half \
                                                             + dot(hat_psi_aux,dot(ANN,hat_psi_aux)))\
                      +H0*g*W_half*delta_h*(h_n1-H) - H0*psi_1_half*X*dWM_half_dt*delta_h.dx(0)))*dx \
                        - (Lw*dWM_half_dt*delta_h*(psi_1_half*I1 + dot(hat_psi_aux,IN)))*ds(1)

    elif dim=="2D":
        a_psi_1 = (H0*delta_h*(Lw-WM_n1)*(psi_1)/(0.5*dt))*dx
        L_psi_1 = -( -H0*delta_h*(Lw-WM_half)*psi_1_half/(0.5*dt) \
                    +(delta_h*(((Lw*Lw)/(2.0*(Lw-WM_half)))*((psi_1_half.dx(0)**2)*M11 +dot(hat_psi_aux.dx(0), \
                                                                                            (2.0*MN1*psi_1_half.dx(0)\
                                                                                             +dot(MNN,hat_psi_aux.dx(0))))))\
                      -((1.0/(Lw-WM_half))*(Lw*Lw)*delta_h.dx(0))*( psi_1_half.dx(0)*(D11*psi_1_half + dot(D1N,hat_psi_aux)) \
                                                                   +dot(hat_psi_aux.dx(0), (DN1*psi_1_half + dot(DNN, hat_psi_aux))))\
                      +(1.0/h_n1)*(delta_h.dx(0)*((1.0/(Lw-WM_half))*h_n1.dx(0)*(Lw*Lw))\
                                   -(delta_h/h_n1)*( (Lw*Lw)*(h_n1.dx(0)**2)/(2.0*(Lw-WM_half))))*(psi_1_half*psi_1_half*S11 \
                                                                                                   + 2.0*dot(hat_psi_aux,SN1)*psi_1_half\
                                                                                                   +dot(hat_psi_aux,dot(SNN,hat_psi_aux)))\
                      -(0.5*delta_h*(Lw-WM_half)*H0*H0/(h_n1**2))*(psi_1_half*psi_1_half*A11 + 2.0*dot(hat_psi_aux,AN1)*psi_1_half \
                                                                   + dot(hat_psi_aux,dot(ANN,hat_psi_aux)))\
                      +H0*g*(Lw-WM_half)*delta_h*(h_n1-H) - H0*psi_1_half*(x_coord-Lw)*dWM_half_dt*delta_h.dx(0)))*dx \
                        - (Lw*dWM_half_dt*delta_h*(psi_1_half*I1 + dot(hat_psi_aux,IN)))*ds(1)


    return a_psi_1, L_psi_1



#----------------------------------------------------------------------------------------------------------------------#
#                                        Step 4 : Update psi_i at time t^{n+1}:                                        #
#______________________________________________________________________________________________________________________#

def WF_hat_psi_SV(b,dim, n_z, Lw, H0, H, WM, x_coord, dt, dWM_dt, dWM_dy, delta_hat_psi, hat_psi, h_n0, psi_1_n0, M11, M1N, MN1, MNN, D11, D1N, DN1, DNN, S11, SN1, SNN, A11, AN1, ANN, I1, IN, B11, B1N, BN1, BNN, C11, C1N, CN1, CNN  ):
    if dim=="3D":
        a_hat_psi =((h_n0/(Lw-WM))*(Lw*Lw+((x_coord-Lw)*dWM_dy)**2)*elem_mult(delta_hat_psi.dx(0), dot(MNN,hat_psi.dx(0)))\
                    +(Lw-WM)*h_n0*elem_mult(delta_hat_psi.dx(1),dot(MNN,hat_psi.dx(1)))\
                    +(x_coord-Lw)*dWM_dy*h_n0*(elem_mult(dot(hat_psi.dx(0),MNN),delta_hat_psi.dx(1))\
                                               +elem_mult(delta_hat_psi.dx(0),dot(MNN,hat_psi.dx(1))))\
                    -( (Lw*Lw+((x_coord-Lw)*(dWM_dy))**2)*h_n0.dx(0)/(Lw-WM) \
                      +(x_coord-Lw)*dWM_dy*h_n0.dx(1))*(elem_mult(delta_hat_psi, dot(DNN.T,hat_psi.dx(0))) \
                                                        + elem_mult(delta_hat_psi.dx(0),dot(DNN,hat_psi))) \
                    -((Lw-WM)*h_n0.dx(1) + (x_coord-Lw)*dWM_dy*h_n0.dx(0))*(elem_mult(delta_hat_psi, dot(DNN.T,hat_psi.dx(1)))\
                                                                            +elem_mult(delta_hat_psi.dx(1),dot(DNN,hat_psi)))\
                    -H0*( (Lw*Lw+((x_coord-Lw)*(dWM_dy))**2)*b.dx(0)/(Lw-WM) \
                      +(x_coord-Lw)*dWM_dy*b.dx(1))*(elem_mult(delta_hat_psi, dot(BNN.T,hat_psi.dx(0))) \
                                                        + elem_mult(delta_hat_psi.dx(0),dot(BNN,hat_psi))) \
                    -H0*((Lw-WM)*b.dx(1) + (x_coord-Lw)*dWM_dy*b.dx(0))*(elem_mult(delta_hat_psi, dot(BNN.T,hat_psi.dx(1)))\
                                                                            +elem_mult(delta_hat_psi.dx(1),dot(BNN,hat_psi)))\
                    +(1.0/h_n0)*( (Lw*Lw+((x_coord-Lw)*(dWM_dy))**2)*(h_n0.dx(0)**2)/(Lw-WM) +(Lw-WM)*(h_n0.dx(1)**2) \
                                 +2.0*h_n0.dx(0)*h_n0.dx(1)*(x_coord-Lw)*dWM_dy)*elem_mult(delta_hat_psi,dot(SNN,hat_psi))\
                    +(2.0*H0/h_n0)*( (Lw*Lw+((x_coord-Lw)*(dWM_dy))**2)*(h_n0.dx(0)*b.dx(0))/(Lw-WM) +(Lw-WM)*(h_n0.dx(1)*b.dx(1)) \
                                 +b.dx(0)*h_n0.dx(1)*(x_coord-Lw)*dWM_dy+b.dx(1)*h_n0.dx(0)*(x_coord-Lw)*dWM_dy)*elem_mult(delta_hat_psi,dot(CNN,hat_psi))\
                    + (H0*H0/h_n0)*((Lw-WM)+(Lw*Lw+((x_coord-Lw)*(dWM_dy))**2)*(b.dx(0)**2)/(Lw-WM) +(Lw-WM)*(b.dx(1)**2) \
                                    +2.0*b.dx(0)*b.dx(1)*(x_coord-Lw)*dWM_dy)*elem_mult(delta_hat_psi,dot(ANN,hat_psi)))
        
        L_hat_psi =-((h_n0/(Lw-WM))*(Lw*Lw+((x_coord-Lw)*dWM_dy)**2)*elem_mult(delta_hat_psi.dx(0), MN1*psi_1_n0.dx(0))\
                     +(Lw-WM)*h_n0*elem_mult(delta_hat_psi.dx(1), MN1*psi_1_n0.dx(1))\
                     +(x_coord-Lw)*dWM_dy*h_n0*(elem_mult(delta_hat_psi.dx(1), MN1*psi_1_n0.dx(0))\
                                                +elem_mult(delta_hat_psi.dx(0),MN1*psi_1_n0.dx(1)))\
                     -( (Lw*Lw+((x_coord-Lw)*(dWM_dy))**2)*h_n0.dx(0)/(Lw-WM) \
                       +(x_coord-Lw)*dWM_dy*h_n0.dx(1))*(elem_mult(delta_hat_psi, D1N*psi_1_n0.dx(0)) \
                                                         + elem_mult(delta_hat_psi.dx(0),DN1*psi_1_n0)) \
                     -((Lw-WM)*h_n0.dx(1) + (x_coord-Lw)*dWM_dy*h_n0.dx(0))*(elem_mult(delta_hat_psi, D1N*psi_1_n0.dx(1))\
                                                                             +elem_mult(delta_hat_psi.dx(1),DN1*psi_1_n0))\
                     -H0*( (Lw*Lw+((x_coord-Lw)*(dWM_dy))**2)*b.dx(0)/(Lw-WM) \
                       +(x_coord-Lw)*dWM_dy*b.dx(1))*(elem_mult(delta_hat_psi, B1N*psi_1_n0.dx(0)) \
                                                         + elem_mult(delta_hat_psi.dx(0),BN1*psi_1_n0)) \
                     -H0*((Lw-WM)*b.dx(1) + (x_coord-Lw)*dWM_dy*b.dx(0))*(elem_mult(delta_hat_psi, B1N*psi_1_n0.dx(1))\
                                                                             +elem_mult(delta_hat_psi.dx(1),BN1*psi_1_n0))\
                     +(1.0/h_n0)*( (Lw*Lw+((x_coord-Lw)*(dWM_dy))**2)*(h_n0.dx(0)**2)/(Lw-WM) +(Lw-WM)*(h_n0.dx(1)**2) \
                                  +2.0*h_n0.dx(0)*h_n0.dx(1)*(x_coord-Lw)*dWM_dy)*elem_mult(delta_hat_psi,SN1*psi_1_n0)\
                     +(2.0*H0/h_n0)*( (Lw*Lw+((x_coord-Lw)*(dWM_dy))**2)*(h_n0.dx(0)*b.dx(0))/(Lw-WM) +(Lw-WM)*(h_n0.dx(1)*b.dx(1)) \
                                  +b.dx(0)*h_n0.dx(1)*(x_coord-Lw)*dWM_dy+h_n0.dx(0)*b.dx(1)*(x_coord-Lw)*dWM_dy)*elem_mult(delta_hat_psi,CN1*psi_1_n0)\
                     + (H0*H0/h_n0)*((Lw-WM)+  (Lw*Lw+((x_coord-Lw)*(dWM_dy))**2)*(b.dx(0)**2)/(Lw-WM) +(Lw-WM)*(b.dx(1)**2) \
                                     +2.0*b.dx(0)*b.dx(1)*(x_coord-Lw)*dWM_dy)*elem_mult(delta_hat_psi,AN1*psi_1_n0))
            
    elif dim=="2D":
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


    hat_psi_BC = -(Lw*dWM_dt*h_n0*elem_mult(delta_hat_psi,IN))
    A_hat = sum((a_hat_psi[ind])*dx for ind in range(0,n_z))
    L_hat = sum((L_hat_psi[ind])*dx for ind in range(0,n_z)) + sum((hat_psi_BC[ind])*ds(1) for ind in range(0,n_z))

    return A_hat, L_hat


