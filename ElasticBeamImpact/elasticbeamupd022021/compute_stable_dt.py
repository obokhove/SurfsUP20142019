# -*- coding: utf-8 -*-
"""
Computes maximum stable time step.
"""

import math
from lib.parameters import parameters as prm

rho = prm['rho']
mu = prm['mu']
lam = prm['lam']
dx = prm['Lx'] / prm['nx']
dz = prm['Lz'] / prm['nz']
dy = prm['Ly'] / prm['ny']

k = 2. * math.pi * math.sqrt( 1./dx**2 + 1./dy**2 + 1./dz**2 )
omega = k * math.sqrt( (mu + 2.*lam)/rho )
dt = 2./omega
print('beam dt = ', dt)