"""
Auxiliary functions.

"""

import firedrake as fd
from . import beam
from . import parameters as prm

import os
import numpy as np

fd.parameters["form_compiler"]["cpp_optimize"] = True

def write_energy_to_file(filename, E_tab, time_tab):
    directory = os.path.dirname(filename)
    if not os.path.exists(directory):
        os.makedirs(directory)
    outfile = open(filename, 'w')
    for ind in range(len(time_tab)):
        outfile.write(str(time_tab[ind]) + '\t' + str(E_tab[ind]) + '\n')
    outfile.close()

def initialize_energies(t_len):
    E = dict()
    E['p'] = np.zeros(t_len)
    E['k'] = np.zeros(t_len)
    E['t'] = np.zeros(t_len)
    return E

def update_energies(B, E, i):
    E['p'][i] = B.E_pot()
    E['k'][i] = B.E_kin()
    E['t'][i] = E['p'][i] + E['k'][i]

def write_energies_to_files(E, t_tab):
    write_energy_to_file("energy/Ep", E['p'], t_tab)
    write_energy_to_file("energy/Ek", E['k'], t_tab)
    write_energy_to_file("energy/Et", E['t'], t_tab)
    
def time_evolution():
    params = prm.dim2nondim(prm.parameters)
    B = beam.Beam(**params)
    E = initialize_energies(len(B.t))
    update_energies(B, E, 0)
    n = 0
    n_modulo = 100
    for i, t in enumerate(B.t):
        if i==0: continue
        B.evolve_time()
        if n%n_modulo == 0:
            print('time = ', t * params['T'])
            B.output_data()
        update_energies(B, E, i)
        n+=1
    write_energies_to_files(E, B.t)
    B.write_raw()