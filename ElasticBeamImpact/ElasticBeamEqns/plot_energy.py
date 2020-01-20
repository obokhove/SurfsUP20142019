# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 19:06:43 2015

@author: mmtjs

Plots energy.
"""

import matplotlib.pyplot as plt
import lib.parameters as prm
import numpy as np

fig = plt.figure()#figsize=(7,5))

params = prm.dim2nondim(prm.parameters)

list_of_files = list()
list_of_files.extend( [("energy/Ep", "Ep"), ("energy/Ek", "Ek"), ("energy/Et", "Et")] )

datalist = [ ( np.loadtxt(filename), label ) for filename, label in list_of_files ]

T = params['T']
J = params['M'] * (params['L'] / T )**2

#print 'T = ', T
#print 'J = ', J

for data, label in datalist:
    plt.plot( T * data[:,0], J * data[:,1], label=label )

plt.legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0.)
plt.title("Energy(time)")
plt.xlabel("time [s]")
plt.ylabel("Energy [J]")
#plt.show()

fig.savefig('figs/E_tot.pdf', bbox_inches="tight")
plt.close(fig)


def plot_relative_energies(t, E_tot1, E_tot2):#, E_tot3):
    rel1 = (E_tot1 - E_tot1[0]) / E_tot1[0]
    rel2 = (E_tot2 - E_tot2[0]) / E_tot2[0]
#    rel3 = (E_tot3 - E_tot3[0]) / E_tot3[0]
#    ratio = rel1 / rel2[::2]
    
    fig = plt.figure()
    plt.xlim( (t[0], t[-1]) )
#    plt.plot(t[::2], ratio)
    plt.plot( t, rel1, ':', color='darkgray', label=r'@ $\Delta t$' )
    plt.plot( t, rel2[:len(rel1)*2:2], 'k', lw=1., label=r'@ $\Delta t / 2 $' )
#    plt.plot( t, rel3[:len(rel1)*4:4], 'k', label=r'@ $\Delta t / 4 $' )
#    plt.loglog( dx_t, L2_E, '--*', label='$L_2$')
#    plt.loglog( dx_t, Linf_E, '--o', label='$L_{\infty}$')
    plt.title(r'$(E(t) - E(0))/E(0)$ as $ \Delta t \rightarrow \Delta t / 2$')
    plt.xlabel('t[s]')
    plt.ylabel(r'$(E(t) - E(0))/E(0)$')
    plt.legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0.)
    fig.savefig('figs/rel_energy_dt_dt2.pdf', bbox_inches='tight')
    plt.close(fig)


#def plot_E_tot(CS, E_tot ):
#    fig = plt.figure()#figsize=(7,5))
#    t = CS.t * CS.parameters['T']
#    J = CS.parameters['M'] * (CS.parameters['L'] / CS.parameters['T'])**2
#    plt.xlim( (t[0], t[-1]) )
#    plt.plot(t, J*E_tot, 'k', lw=1.5, label='total')
#    plt.plot(t, J*CS.W.E_tot, '', color='darkgray', lw=2.5, label='tot water') # color='sage'
#    plt.plot(t, J*CS.W.E_pot, ':', color='darkgray', lw=2.5, label='pot water')
#    plt.plot(t, J*CS.W.E_kin, '--', color='darkgray', lw=2.5, label='kin water')
#    plt.plot(t, J*CS.B.E_tot, 'k', label='tot beam')
#    plt.plot(t, J*CS.B.E_pot, 'k:',  label='pot beam')
#    plt.plot(t, J*CS.B.E_kin, 'k--', label='kin beam')
#    if CS.W.flags['Wexact']:
#        plt.plot(t, J*CS.W.E_kin_exact, 'r', label='W exact kin')
#        plt.plot(t, J*CS.W.E_pot_exact, '--r', label='W exact pot')
#        plt.plot(t, J*CS.W.E_tot_exact, ':r', label='W exact tot')
#    plt.title('Energy E(t)')
#    plt.xlabel('t [s]'); plt.ylabel('E [J]')
#    plt.legend(bbox_to_anchor=(1.005, 1), loc='upper left', borderaxespad=0.)
#    fig.savefig('figs/E_tot.pdf', bbox_inches="tight")
#    plt.close(fig)