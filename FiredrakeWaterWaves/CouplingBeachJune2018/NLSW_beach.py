# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 14:28:55 2015

@author: mmfg
"""


""" 
    Coupling : On average h et hu dans la derniere cellule de PF, et on le met dans la premiere cellule de FV. Ensuite on prend le flux at x0+1/2 de FV pour l'utiliser a la boundary de PF (car besoin d'un flux). Donc PF et SW ont une cellule en commun, le changer dans le rapport (certainement ce que disait Onno).
    """
import numpy as np
from firedrake import *

def U_half(h, U):
    u = 0*U[0,:]
    [Ind] = np.where(U[0,:]!=0)
    u[Ind]=U[1,Ind]/U[0,Ind]

    u_half = 0*np.eye(2,len(U[0,:]))
    u_half[0,:]=h
    u_half[1,:]=h*u
    return [u,u_half]


def right_wave_speed(ul, ur, hl, hr, g):
    Sr = np.maximum(ul+ np.sqrt(g*hl), ur+ np.sqrt(g*hr))
    return Sr

def left_wave_speed(ul, ur, hl, hr, g):
    Sl = np.minimum(ul- np.sqrt(g*hl), ur- np.sqrt(g*hr))
    return Sl

def flux_half(U_half,g):
    F=0*np.eye(2,len(U_half[0,:]))
    F[0,:]=U_half[1,:]
    [Ind] = np.where(U_half[0,:]!=0)
    F[1,Ind] = U_half[1,Ind]*U_half[1,Ind]/U_half[0,Ind] + 0.5*g*U_half[0,Ind]**2
    return F


def HLL_flux(Fl, Fr, Sl, Sr, Ul, Ur):
    F = 0*np.eye(2,len(Sl))
    [Ind] = np.where(Sl>0)
    F[:,Ind] = Fl[:,Ind]        #Fl
    [Ind] = np.where((Sl<=0) & (0<=Sr) & (Sl!=Sr))
    F[0,Ind] = (Sr[Ind]*Fl[0,Ind]-Sl[Ind]*Fr[0,Ind] + Sl[Ind]*Sr[Ind]*(Ur[0,Ind]-Ul[0,Ind]))/(Sr[Ind]-Sl[Ind])
    F[1,Ind] = (Sr[Ind]*Fl[1,Ind]-Sl[Ind]*Fr[1,Ind] + Sl[Ind]*Sr[Ind]*(Ur[1,Ind]-Ul[1,Ind]))/(Sr[Ind]-Sl[Ind])
    [Ind] = np.where(Sr<0)
    F[:,Ind] = Fr[:,Ind]
    return F


def init_FV_solutions(N, H0, bk, i_half, h_bc, hu_bc, h, u):
    # --- Initial values
    Nvol=N
    ul=0.0
    ur=0.0
    Al=H0-bk[1:i_half+1]          # indices start from 0 (interval ouvert)
    Ar=H0-bk[i_half+1:Nvol+1]    # last index = Nvol+1 (but not taken)
    Al[Al<0.0]=0.0
    Ar[Ar<0.0]=0.0

    U=0*np.eye(2,Nvol+2)
    # h
    U[0,1:i_half+1]=np.copy(Al)
    U[0,i_half+1:Nvol+1]=np.copy(Ar)      # indices go from 0 to Nvol+1 to give Nvol+2 values
    U[0,0]=h_bc                             # hr=hl
    U[0,Nvol+1]=U[0,Nvol]
    # hu
    U[1,1:i_half+1]=Al*ul
    U[1,i_half+1:Nvol+1]=Ar*ur   # indices go from 0 to Nvol+1 to give Nvol+2 values
    U[1,0]= hu_bc
    U[1,Nvol+1]=-U[1,Nvol]    # hrur =-hl ul : reflection on the right wall
    
    return U
    
def FV_solutions(U, N):
    Nvol=N
    u_sol = 0*np.eye(1,Nvol+2)
    h_sol = 0*np.eye(1,Nvol+2)
    u_sol[0,np.where(U[0,:]!=0)]= U[1,np.where(U[0,:]!=0)]/U[0,np.where(U[0,:]!=0)]
    h_sol = U[0,:]
    hu_sol = U[1,:]
    return h_sol, u_sol, hu_sol
    
def solve_FV(Nvol, d_x, dt, sw_beach, g, h_bc, hu_bc, h_fv, hu_fv, hu_fe):
    g_tilde=g
    
    #------------ Topography -----------#
    bk = 0*np.eye(1,Nvol+1)
    bk[0,1:Nvol+1] = sw_beach.dat.data[:]
    bk[0,0] = bk[0,1]-(bk[0,2]-bk[0,1])

    #---------- Solutions U^n: ----------#
    U = 0*np.eye(2,Nvol+1)
    U[0,1:Nvol+1] = h_fv.dat.data[0:Nvol]
    U[1,1:Nvol+1] = hu_fv.dat.data[0:Nvol]
    U[1,Nvol]=U[1,Nvol-1]
    U[0,Nvol]=U[0,Nvol-1]
    
    #-------- Boundary conditions -------#
    U[0,0]=h_bc
    U[1,0]=hu_bc
    
    #-------- Initialise U^{n+1} --------#
    U_next=np.copy(U)
    
    
    ############################################################################
    #                                 Update h                                 #
    ############################################################################
    k = range(1,Nvol-1)
    k_plus = range(2,Nvol)
    k_minus = range(0,Nvol-2)
    Uk = np.copy(U[:,k])             # Uk
    Uk_plus = np.copy(U[:,k_plus])   # U_{k+1}
    Uk_minus = np.copy(U[:,k_minus]) # U_{k-1}
    
    #--------------------------- Non-negative depth ---------------------------#
    bk_half_r = np.maximum(bk[0,k],bk[0,k_plus])                     # b_{k+1/2}
    bk_half_l = np.maximum(bk[0,k],bk[0,k_minus])                    # b_{k-1/2}
    h_plus_r = np.maximum(Uk_plus[0,:]+bk[0,k_plus]-bk_half_r,0)   # h_{k+1/2^+}
    h_plus_l = np.maximum(Uk[0,:]+bk[0,k]-bk_half_r,0)             # h_{k+1/2^-}
    h_minus_r = np.maximum(Uk[0,:]+bk[0,k]-bk_half_l,0)            # h_{k-1/2^+}
    h_minus_l = np.maximum(Uk_minus[0,:]+bk[0,k_minus]-bk_half_l,0)# h_{k-1/2^-}

    #----- Left and right values of U at the interface k+/-1/2 -----#
    [u_plus_l,U_plus_l] = U_half(h_plus_l, Uk)            # U(k+1/2)-
    [u_plus_r,U_plus_r] = U_half(h_plus_r, Uk_plus)       # U(k+1/2)+
    [u_minus_l,U_minus_l] = U_half(h_minus_l, Uk_minus)   # U(k-1/2)-
    [u_minus_r,U_minus_r] = U_half(h_minus_r, Uk)         # U(k-1/2)+
    
    # Left and right values of F at the interface k+/-1/2 #
    Fl_minus =flux_half(U_minus_l,g_tilde)  # Fl(U_{k-1/2})
    Fr_minus = flux_half(U_minus_r,g_tilde) # Fr(U_{k-1/2})
    Fl_plus = flux_half(U_plus_l,g_tilde)   # Fl(U_{k+1/2})
    Fr_plus = flux_half(U_plus_r,g_tilde)   # Fr(U_{k+1/2})
    
    #--------------------- Left and right speeds at each interface of cell k ----------------------#
    Sl_minus = left_wave_speed(u_minus_l,u_minus_r,h_minus_l,h_minus_r,g_tilde) # Ul/r = U(k-1/2)-/+
    Sr_minus = right_wave_speed(u_minus_l,u_minus_r,h_minus_l,h_minus_r,g_tilde)# Ul/r = U(k-1/2)-/+
    Sl_plus = left_wave_speed(u_plus_l,u_plus_r,h_plus_l,h_plus_r,g_tilde)      # Ul/r = U(k+1/2)-/+
    Sr_plus = right_wave_speed(u_plus_l,u_plus_r,h_plus_l,h_plus_r,g_tilde)     # Ul/r = U(k+1/2)-/+
    
    #---------------------------------- HLL fluxes ----------------------------------#
    # F(k-1/2) = F(Fl,Fr) ; Fl = F(U(k-1/2)-) ; Fr = F(U(k-1/2)+)
    F_minus = HLL_flux(Fl_minus, Fr_minus,  Sl_minus, Sr_minus, U_minus_l, U_minus_r )
    # F(k+1/2) = F(Fl,Fr) ; Fl = F(U(k+1/2)-) ; Fr = F(U(k+1/2)+)
    F_plus  = HLL_flux(Fl_plus, Fr_plus, Sl_plus, Sr_plus, U_plus_l, U_plus_r)
    
    #--------------------- Topography ---------------------#
    Sk = 0*np.eye(2,len(k))
    
    #-------------------------------- Update h --------------------------------#
    U_next[0,k] = Uk[0,:] - dt*(F_plus[0,:] - F_minus[0,:])/d_x + dt*Sk[0,:]/d_x
    
    # Boundary condition
    U_next[0,-1]=U_next[0,-2] # hr = hr



    ############################################################################
    #                                 Update hu                                #
    ############################################################################
    #--- Left and Right values of U^* ---#
    Hk = np.copy(U_next[0,k])
    Hk_plus =np.copy(U_next[0,k_plus])
    Hk_minus = np.copy(U_next[0,k_minus])
        
    #--------------------- Non-negative depth ---------------------#
    h_next_plus_r = np.maximum(Hk_plus+bk[0,k_plus]-bk_half_r,0)
    h_next_plus_l = np.maximum(Hk+bk[0,k]-bk_half_r,0)
    h_next_minus_r = np.maximum(Hk+bk[0,k]-bk_half_l,0)
    h_next_minus_l = np.maximum(Hk_minus+bk[0,k_minus]-bk_half_l,0)



#
    #-- Left and right values of U at the interface k+1/2 --#
    U_plus_l[0,:] = h_next_plus_l         # h_{k+1/2^-}^{n+1}
    U_plus_r[0,:] = h_next_plus_r         # h_{k+1/2^+}^{n+1}
    U_minus_l[0,:] = h_next_minus_l       # h_{k-1/2^-}^{n+1}
    U_minus_r[0,:] = h_next_minus_r       # h_{k-1/2^+}^{n+1}
    
    #- Left and right values of F at the interface k+/-1/2 -#
    Fl_minus =flux_half(U_minus_l,g_tilde)    # Fl(U_{k-1/2})
    Fr_minus = flux_half(U_minus_r,g_tilde)   # Fr(U_{k-1/2})
    Fl_plus = flux_half(U_plus_l,g_tilde)     # Fl(U_{k+1/2})
    Fr_plus = flux_half(U_plus_r,g_tilde)     # Fr(U_{k+1/2})
        
    #---------------------------- Left and right speeds at each interface of cell k -----------------------------#
    Sl_minus = left_wave_speed(u_minus_l, u_minus_r, h_next_minus_l, h_next_minus_r,g_tilde)  # Ul/r = U(k-1/2)-/+
    Sr_minus = right_wave_speed(u_minus_l, u_minus_r, h_next_minus_l,h_next_minus_r, g_tilde) # Ul/r = U(k-1/2)-/+
    Sl_plus = left_wave_speed(u_plus_l, u_plus_r, h_next_plus_l, h_next_plus_r, g_tilde)      # Ul/r = U(k+1/2)-/+
    Sr_plus = right_wave_speed(u_plus_l, u_plus_r, h_next_plus_l, h_next_plus_r, g_tilde)     # Ul/r = U(k+1/2)-/+


    #---------------------------------- HLL fluxes ----------------------------------#
    # F(k-1/2) = F(Fl,Fr) ; Fl = F(U(k-1/2)-) ; Fr = F(U(k-1/2)+)
    F_minus = HLL_flux(Fl_minus, Fr_minus,  Sl_minus, Sr_minus, U_minus_l, U_minus_r )
    # F(k+1/2) = F(Fl,Fr) ; Fl = F(U(k+1/2)-) ; Fr = F(U(k+1/2)+)
    F_plus  = HLL_flux(Fl_plus, Fr_plus, Sl_plus, Sr_plus, U_plus_l, U_plus_r)

    #--------------------- Topography ---------------------#
    Sk[1,:]=0.5*g_tilde*h_next_plus_l**2-0.5*g_tilde*h_next_minus_r**2
    
    #-------------------------------- Update hu -------------------------------#
    U_next[1,k] = Uk[1,:] - dt*(F_plus[1,:] - F_minus[1,:])/d_x + dt*Sk[1,:]/d_x

    # BC for DW: the flux at the left boundary
    hu_fe.vector().set_local(F_minus[0,0])

    #--------- Boundary conditions --------#
    U_next[1,-1]=-U_next[1,-2] # -hur = -hul
    
    #---- Update -----#
    U = np.copy(U_next)
    
    #-------- Update the solutions -------#
    h_fv.vector().set_local(U[0,1:Nvol+1])
    hu_fv.vector().set_local(U[1,1:Nvol+1])
    
    #---- Compute the energy ----#
    hu_square = 0.0*np.eye(1,Nvol)
    h_square = 0.0*np.eye(1,Nvol)
    kk = range(1,Nvol)
    kk_minus = range(0,Nvol-1)
    huu = 0.0*np.eye(1,len(U[0,:]))
    [Ind] = np.where(U[0,:]>=1e-9)
    huu[0,Ind] =U[1,Ind]*U[1,Ind]/U[0,Ind]
    hu_square[0,kk_minus]=huu[0,kk]
    h_square[0,kk_minus]=U[1,kk]*U[1,kk]
    E_sw = d_x*(0.5*sum(hu_square[0,:]) + 0.5*g*sum(h_square[0,:]))
    return h_fv, hu_fv, U, hu_fe, E_sw


