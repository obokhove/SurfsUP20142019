#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 12:17:37 2017

@author: mmtjs

Program solving compressible van der Waals like fluid equations.
"""

import timeit
import sys

from firedrake import *
import numpy as np

start = timeit.default_timer()

# parameters
c_w = 120 # m/s
c_a = 100 # m/s
c_m = 2. # m/s
H0 = 40. # m
H = 80. # m
Lx = 80. # m
p_a = 1e5 # Pa
rho_w = 1000 # kg/m3
rho_a = p_a / c_a**2 # kg/m3
g = 9.81 # m/s2
k = pi/Lx # first mode of standing wave
omega = 0.600739
dt = 1e-3 # s
t_end = 100. # s
two_phase = False # 2- or 3-phase model

rho_s = (p_a - c_m**2 * (rho_a+rho_w)/2 )/( c_a**2-c_m**2 )
rho_ss = (c_w**2 * rho_w - c_m**2 * (rho_a+rho_w)/2 )/( c_w**2-c_m**2 )
Hs = H0 + c_m**2/g * ln(rho_ss/rho_s)
print('Hs - H0 = ', Hs-H0)
print('rho_a = ', rho_a)
print('rho_w = ', rho_w)
print('rho_s = ', rho_s)
print('rho_ss = ', rho_ss)

fps_process = 20
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
Nx = 20
Nz = 20
prog = 1.1 # progression
d_x = Lx / Nx # mind nameing conflict with UFL dx
dz_min = H0 * (1.-prog)/(1.-prog**Nz)
print( "dx = ", L * d_x )
print( "min dz = ", L * dz_min ) 
st_dt = dz_min/pi/c_w
print( "stable dt SV = ", T * st_dt ) # for SV/SE
k_max = pi / dz_min
dt_max = 2. * pi / sqrt( k_max * tanh(k_max * H0) )
print( "max dt water wave = ", T * dt_max )
T_standing_wave = 2.*pi/sqrt(k*tanh(k*H0))
print( 'T standing wave = ', T * T_standing_wave ) # nondim

#quit()

#mesh = UnitSquareMesh(Nx, Nz)
mesh = Mesh("mesh_coarse.msh")
V = FunctionSpace(mesh, "CG", 1)
rho = Function(V, name="rho")
rho_0 = Function(V, name="rho_h")
phi = Function(V, name="phi")

v = TestFunction(V)
phi_next = Function(V)
rho_next = Function(V)

trial = TrialFunction(V)

x, z = SpatialCoordinate(mesh)

ampl = 0.1 * H0
eta0 = ampl * cos(k*x)

rho_lin = Function(V, name="rho_lin")
rho_lin.interpolate( conditional(z>(eta0+H0), 0., rho_w) )

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

def mi(a, b):
    return conditional(a<b, a, b)

def ma(a, b):
    return conditional(a>b, a, b)

def Qsharp(rho):
    return c_a**2 * ln( conditional(rho<rho_a, rho, rho_a)
                       / conditional(rho_0<rho_a, rho_0, rho_a) ) \
        + c_w**2 * ln( conditional( rho>rho_w, rho, rho_w )
                      / conditional( rho_0>rho_w, rho_0, rho_w ) )

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
#solver_parameters = {
#                    'ksp_atol': 1e-30,
#                   'ksp_rtol': 1e-9,
#                   'ksp_divtol': 1e4
#                   }

solver_parameters={
    'mat_type':'aij',
    'ksp_type': 'preonly',
    'pc_type': 'lu', 
    'pc_factor_mat_solver_package': 'mumps',
    'ksp_rtol': 1e-14
}

a_rho = trial * (v - dt * inner( grad(v), grad(phi) ) )*dx
L_rho = v*rho*dx
LVP_rho = LinearVariationalProblem( a_rho, L_rho, rho )
LVS_rho = LinearVariationalSolver( LVP_rho, solver_parameters=solver_parameters )

a_phi = trial*v*dx
L_phi = -( (-phi) + dt/2.*inner(grad(phi), grad(phi)) + dt*Q(rho) )*v*dx
LVP_phi = LinearVariationalProblem( a_phi, L_phi, phi )
LVS_phi = LinearVariationalSolver( LVP_phi, solver_parameters=solver_parameters )

outfile_phi = File("results/phi.pvd")
outfile_rho = File("results/rho.pvd")

outfile_phi.write(phi)
outfile_rho.write(rho, rho_lin, time=0.)

t=0.
n=0

while t<=t_end:
    LVS_rho.solve()
    LVS_phi.solve()
    t += dt
    n += 1
    if n%n_modulo == 0:
        print('t = ', T*t)
        eta = eta0 * cos(omega*t)
        rho_lin.interpolate( conditional(z>(eta+H0), 0., rho_w) )
        outfile_phi.write(phi, time=T*t)
        outfile_rho.write(rho, rho_lin, time=T*t)
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