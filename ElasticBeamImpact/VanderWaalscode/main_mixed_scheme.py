#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 10:45:46 2019

@author: mmtjs

Program solving compressible gas equations.
Mixed space, between modified midpoint and Euler scheme.
"""

import timeit

from firedrake import *
import numpy as np

start = timeit.default_timer()

# parameters
c_w = 120 # m/s
c_a = 100 # m/s
c_m = 3. # m/s
H0 = 40. # m
H = 80. # m
Lx = 80. # m
p_a = 1e5 # Pa
rho_w = 1000 # kg/m3
rho_a = p_a / c_a**2 # kg/m3
g = 9.81 # m/s2
k = pi/Lx # first mode of standing wave
omega = 0.600739
dt = 1e-2 # s
t_end = 100. # s
two_phase = True # 2- or 3-phase model
a = 0.55 # scheme impliciteness parameter

rho_s = (p_a - c_m**2 * (rho_a+rho_w)/2 )/( c_a**2-c_m**2 )
rho_ss = (c_w**2 * rho_w - c_m**2 * (rho_a+rho_w)/2 )/( c_w**2-c_m**2 )
Hs = H0 + c_m**2/g * ln(rho_ss/rho_s)
print('Hs - H0 = ', Hs-H0)
print('rho_a = ', rho_a)
print('rho_w = ', rho_w)
print('rho_s = ', rho_s)
print('rho_ss = ', rho_ss)

fps_process = 20 # save this many times per physical second of simulated process
n_modulo = int(round(1/dt,0))//fps_process

# nondim units
L = H
T = sqrt(L/g)

dt /= T
t_end /= T
c_w /= L/T
c_a /= L/T
c_m /= L/T
rho_scale = rho_w
rho_a /= rho_scale
rho_w /= rho_scale
rho_s /= rho_scale
rho_ss /= rho_scale
H0 /= L
Hs /= L
H /= L
Lx /= L
omega *= T
k *= L


# parameters from mesh.geo
Nx = 40
Nz = 80
#Nx = 20 # for coarse
#Nz = 20
prog = 1.07 # progression
d_x = Lx / Nx # mind nameing conflict with UFL dx
dz_min = H0 * (1.-prog)/(1.-prog**Nz)
print( "dx = ", L * d_x )
print( "min dz = ", L * dz_min ) 
st_dt = dz_min/pi/c_w
print( "stable dt SE = ", T * st_dt ) # for SV/SE
k_max = pi / dz_min
dt_max = 2. * pi / sqrt( k_max * tanh(k_max * H0) )
print( "max dt water wave = ", T * dt_max )
T_standing_wave = 2.*pi/sqrt(k*tanh(k*H0))
print( 'T standing wave = ', T * T_standing_wave ) # nondim


#mesh = UnitSquareMesh(Nx, Nz)
#mesh = Mesh("mesh.msh")
mesh = Mesh("mesh_coarse.msh")
#mesh = Mesh("mesh_optimal.msh")
V = FunctionSpace(mesh, "CG", 1)
W = V*V

var_rho, var_phi = TestFunctions(W)
rho_0 = Function(V, name="rho_0") # hydrostatic state
#phihat = Function(V, name="phi_hat")
#phihatoverc2 = Function(V)
phi_lin = Function(V, name="phi_lin")
rho_lin = Function(V, name="rho_lin")

w = Function(W)
w_next = Function(W)
phi, rho = w.split() # mind which type of "split"!!!
phi_next, rho_next = split(w_next) # mind which type of "split"!!!
phi.rename("phi")
rho.rename("rho")


x, z = SpatialCoordinate(mesh)

ampl = 0.05 * H0
eta0 = ampl * cos(k*x)

if two_phase:
    rho_0.interpolate( conditional(z>H0, rho_a * exp(-(z-H0)/c_a**2), rho_w * exp(-(z-H0)/c_w**2) ) )
    rho.interpolate( conditional(z>H0+eta0, rho_a * exp(-(z-H0-eta0)/c_a**2),
                                         rho_w * exp(-(z-H0-eta0)/c_w**2) ) )
else:
    rho_0.interpolate( conditional(z>H0,
                       conditional(z>=Hs,
                        rho_s * exp(-(z-Hs)/c_a**2),
                        rho_ss * exp(-(z-H0)/c_m**2) ),
                        rho_ss * exp(-(z-H0)/c_w**2) ) )
    rho.interpolate( conditional(z>H0+eta0,
                       conditional(z>=Hs+eta0,
                        rho_s * exp(-(z-Hs-eta0)/c_a**2),
                        rho_ss * exp(-(z-H0-eta0)/c_m**2) ),
                        rho_ss * exp(-(z-H0-eta0)/c_w**2) ) )

#rho.assign(rho_0) # just for testing the hydrostatic case

def mi(a, b): # select smaller
    return conditional(a<b, a, b)

def ma(a, b): # select larger
    return conditional(a>b, a, b)

def sig(a): # sign function
    return conditional(a>0, 1, conditional(a<0,-1,0) )

def Qsharp(rho):
    return c_a**2 * ln( mi(rho, rho_a) / mi(rho_0, rho_a) ) \
         + c_w**2 * ln( ma(rho, rho_w) / ma(rho_0, rho_w) )

def Qsmooth(rho):
    tmp1 = c_a**2 * ln( mi(rho, rho_s ) / mi(rho_0, rho_s ) ) \
          +c_w**2 * ln( ma(rho, rho_ss) / ma(rho_0, rho_ss) )
    tmp2 = c_m**2 * conditional( rho_0<=rho,
                      conditional(And(rho_0<rho_ss, rho_s<rho), 
                                  ln( mi(rho,rho_ss) / ma(rho_0,rho_s) ), 0.),
                      conditional( And( rho<rho_ss, rho_s<rho_0 ),
                                  ln( ma(rho,rho_s) / mi(rho_0,rho_ss) ), 0. ) )
    return tmp1 + tmp2

def Q(rho):
    if two_phase:
        return Qsharp(rho)
    else:
        return Qsmooth(rho)


# equations
F_phi = ( (phi_next - phi)/dt \
         +     a*( 1./2.*inner(grad(phi_next), grad(phi_next)) + Q(rho_next) ) \
         + (1.-a)*( 1./2.*inner(grad(phi), grad(phi)) + Q(rho) ) ) * var_rho * dx
F_rho = ( var_phi*(rho_next-rho)/dt \
         - a*( rho_next*inner( grad(var_phi), grad(phi_next) ) ) \
         - (1.-a)*( rho*inner( grad(var_phi), grad(phi) ) ) ) * dx
F = F_phi + F_rho


outfile_phi = File("./results/phi.pvd")
outfile_rho = File("./results/rho.pvd")

t=0.
n=0

#phi_lin.interpolate( L**2/T * phihat*cos(k*x)*cos(omega*t) )
#rho_lin.interpolate( rho_scale * rho_0*(1. + omega*phihatoverc2*cos(k*x)*sin(omega*t)) )

phi_dim = Function(V, name="phi")
rho_dim = Function(V, name="rho")
phi_dim.interpolate( L**2/T * phi )
rho_dim.interpolate( rho_scale * rho )
outfile_phi.write( phi_dim, time = T * t )
outfile_rho.write( rho_dim, time = T * t )
#outfile_phi.write( phi_dim, phi_lin, time = T * t )
#outfile_rho.write( rho_dim, rho_lin, time = T * t )
#outfile_phi.write( phi_lin, time = T * t )
#outfile_rho.write( rho_lin, time = T * t )


#solver_parameters={
#        'snes_type': 'newtonls',
##        'snes_type': 'nrichardson',
##        'snes_type': 'test',
##        'ksp_type': 'gmres', # with a restart (ksp_gmres_restart) of 30
##        'snes_rtol': 1e-8,
##        'snes_atol': 1e-50,
##        'snes_stol': 1e-8,
##        'snes_max_it': 100,
#        'ksp_rtol': 1e-7,
##        'ksp_atol': 1e-50,
##        'ksp_divtol': 1e4,
##        'ksp_max_it': 10000,
##        'pc_type': 'lu' # (Jacobi preconditioning for mixed problems)
#        'mat_type':'aij',
#}

solver_parameters = {
  "snes_type": "newtonls",
#  "snes_type": "nrichardson",
  "ksp_type": "preonly",
  "mat_type": "aij",
  "pc_type": "lu",
  "pc_factor_mat_solver_type": "mumps",
#  "snes_monitor": True,
#  "snes_linesearch_monitor": True,
#  "snes_converged_reason": True,
#  "ksp_converged_reason": True
}

problem = NonlinearVariationalProblem(F, w_next)
solver = NonlinearVariationalSolver(problem, solver_parameters=solver_parameters)

while t <= t_end:
    solver.solve()
#    solve(F == 0, w_next, solver_parameters=solver_parameters)
    w.assign(w_next)
    t += dt
    n += 1
    if n%n_modulo == 0:
        print('t = ', T * t)
        phi_dim.interpolate( L**2/T * phi )
        rho_dim.interpolate( rho_scale * rho )
        outfile_phi.write( phi_dim, time = T * t )
        outfile_rho.write( rho_dim, time = T * t )
#        phi_lin.interpolate( L**2/T * phihat*cos(k*x)*cos(omega*t) )
#        rho_lin.interpolate( rho_scale * rho_0*(1. + omega*phihatoverc2*cos(k*x)*sin(omega*t)) )
#        outfile_phi.write( phi_lin, time = T * t )
#        outfile_rho.write( rho_lin, time = T * t )
#        outfile_phi.write( phi_dim, phi_lin, time = T * t )
#        outfile_rho.write( rho_dim, rho_lin, time = T * t )
        if np.isnan(rho.at(0,0)):
            print('NaN spotted!!!')
            sys.exit()


stop = timeit.default_timer()
total_run_time_sec = stop - start
print('')
print('Summary:')
print('T wave = ', T * 2.*pi/sqrt(k*tanh(k*H0)) )
print('dt = ', T * dt)
print('t_end = ', T * t_end)
print('runtime = ', int(total_run_time_sec)//60//60, ' h ',
      int(total_run_time_sec)//60, ' min ', int(total_run_time_sec)%60, ' sec')