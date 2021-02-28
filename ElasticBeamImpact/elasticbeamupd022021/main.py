"""
Created on Thu Oct 22 18:54:28 2015

@author: mmtjs

Main file to run a simulation of the nonlinear elastic beam in FD.
"""

import timeit
from lib import functions as fn

start = timeit.default_timer()

fn.time_evolution()

stop = timeit.default_timer()
total_run_time_sec = stop - start
print('runtime = ', int(total_run_time_sec)//60, ' mins ', int(total_run_time_sec)%60, ' sec')