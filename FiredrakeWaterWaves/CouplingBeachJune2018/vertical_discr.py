"""
    Created on Thu May 12 2016
    
    @author: Floriane Gidel
    """

import numpy as np
from sympy import *


def varphi_expr(i, n, H0):
    z=Symbol('z')
    k = Symbol('k')
    sigma_k = (n-k)*H0/n
    sigma = lambdify(k,sigma_k,"numpy")
    varphi_z = (Product((z-sigma(k))/(sigma(i)-sigma(k)),(k, 0,i-1))*Product((z-sigma(k))/(sigma(i)-sigma(k)),(k, i+1,n))).doit() #this avoids to evaluate the case k=i.
    return varphi_z

def deriv_varphi_expr(varphi_expr):
    z = Symbol('z')
    deriv = diff(varphi_expr,z)
    return deriv
def WM_i(i,n,H0):
    z=Symbol('z')
    expr_WM = varphi_expr(i,n,H0)
    WM = integrate(expr_WM, (z,0,H0))
    return WM

def M_ij(i,j,n,H0):
    z=Symbol('z')
    expr_M = varphi_expr(i,n,H0)*varphi_expr(j,n,H0)
    M = integrate(expr_M, (z,0,H0))
    return M

def A_ij(i,j,n,H0):
    z=Symbol('z')
    expr_A = deriv_varphi_expr(varphi_expr(i,n,H0))*deriv_varphi_expr(varphi_expr(j,n,H0))
    A = integrate(expr_A, (z,0,H0))
    return A

def D_ij(i,j,n,H0): # ATTENTION NON SYMMETRIC !
    z=Symbol('z')
    expr_D = z*varphi_expr(i,n,H0)*deriv_varphi_expr(varphi_expr(j,n,H0))
    D = integrate(expr_D, (z,0,H0))
    return D

def G_i(i,n,H0):
    z=Symbol('z')
    expr_G = z*deriv_varphi_expr(varphi_expr(i,n,H0))
    G = integrate(expr_G, (z,0,H0))
    return G

def S_ij(i,j,n,H0):
    z=Symbol('z')
    expr_D = z*z*deriv_varphi_expr(varphi_expr(i,n,H0))*deriv_varphi_expr(varphi_expr(j,n,H0))
    D = integrate(expr_D, (z,0,H0))
    return D

def I_i(i,n,H0):
    z=Symbol('z')
    expr_I = varphi_expr(i,n,H0)
    I = integrate(expr_I, (z,0,H0))
    return I

