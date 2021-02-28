"""
Created on Tue Mar 17 16:12:11 2015

@author: mmtjs

Module with parameters for the beam wave interaction problem.
"""

#import firedrake as fd
import math

parameters = { # SI units
't_end': 1.8, # s
'dt': 0.00016,
'g': 9.8, # gravitational acceleration
'nx': 4,
'ny': 4,
'nz': 20,
'Lx': 2.,
'Ly': 2.,
'Lz': 20.,
# steel parameters
'rho': 7700., # kg/m^3
'lam': 1e7, #9.695e10, # N/m^2  - first Lame constant 1e7 before
'mu': 1e7, #7.617e10, # N/m^2  - second Lame constant 1e7 before
'nonlin': True # if nonlinear model
}

def dim2nondim(parameters):
    L = parameters['Lz']
    g = parameters['g']
    rho = parameters['rho']
    U = math.sqrt(g*L)
    T = L/U
    M = rho*L**3
    nondim_parameters = parameters.copy()
    nondim_parameters['L'] = L
    nondim_parameters['T'] = T
    nondim_parameters['U'] = U
    nondim_parameters['M'] = M
    nondim_parameters['t_end'] /= T
    nondim_parameters['dt'] /= T
    nondim_parameters['rho'] /= rho
    nondim_parameters['Lx'] /= L
    nondim_parameters['Ly'] /= L
    nondim_parameters['Lz'] /= L
    nondim_parameters['lam'] /= g*rho*L
    nondim_parameters['mu'] /= g*rho*L
    nondim_parameters['T'] = T
    return nondim_parameters