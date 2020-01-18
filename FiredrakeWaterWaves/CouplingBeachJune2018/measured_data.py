
import os.path
from firedrake import *

def run_number(test_case):
    if test_case == "111":
        run_nr = '124'
    elif test_case == "112":
        run_nr = '126'
    elif test_case == "113":
        run_nr = '128'
    elif test_case == "121":
        run_nr = '130'
    elif test_case == "122":
        run_nr = '132'
    elif test_case == "123":
        run_nr = '134'
    elif test_case == "131":
        run_nr = '136'
    elif test_case == "132":
        run_nr = '138'
    elif test_case == "21":
        run_nr = '140'
    elif test_case == "221":
        run_nr = '142'
    elif test_case == "222":
        run_nr = '144'
    elif test_case == "223":
        run_nr = '146'
    elif test_case == "23":
        run_nr = '148'
    return run_nr


def load_wavemaker(measurement_path):
    wm_motion = open(os.path.join(measurement_path, 'PistonMotion.dat'))
    lst=[]
    for line in wm_motion:
        lst+=[line.split()]
    wm_data = [float(x[1]) for x in lst]                   # measured motion

    wm_velocity = open(os.path.join(measurement_path, 'PistonVelocity.dat'))
    lst = []
    for line in wm_velocity:
        lst+=[line.split()]
    t_data = [float(x[0]) for x in lst]                      # measured time
    wm_vel_data = [float(x[1]) for x in lst]                      # velocity
    return wm_data, wm_vel_data, t_data


def interpolate_wavemaker(wm_data, wm_vel_data, t_data, t, dt, Lw):

    WM_expr = Expression("((wm2*(t-t1) - wm1*(t-t2))/(t2-t1))*0.5*(1.0+copysign(1.0,Lw-x[0]))",\
                         wm2=wm_data[1], wm1=wm_data[0], t1=t_data[0], t2=t_data[1], t=t, Lw=Lw)

    dWM_expr = Expression("((dwm2*(t-t1) - dwm1*(t-t2))/(t2-t1))*0.5*(1.0+copysign(1.0,Lw-x[0]))",\
                          dwm2=wm_vel_data[1], dwm1=wm_vel_data[0], t1=t_data[0], t2=t_data[1], t=t, Lw=Lw)

    return WM_expr, dWM_expr

def probe_location(res_dw):
    x1 = 15.002
    x2 = 17.086
    x3 = 19.040
    x4 = 20.015
    x5 = 21.084
    x6 = 22.022
    x7 = 23.159

    Ind_1 = int(x1/res_dw)
    Ind_2 = int(x2/res_dw)
    Ind_3 = int(x3/res_dw)
    Ind_4 = int(x4/res_dw)
    Ind_5 = int(x5/res_dw)
    Ind_6 = int(x6/res_dw)
    Ind_7 = int(x7/res_dw)

    return Ind_1, Ind_2, Ind_3, Ind_4, Ind_5, Ind_6, Ind_7

def probe_files(save_path):
    x1_file = open(os.path.join(save_path, 'probe1.txt'), 'w')
    x2_file = open(os.path.join(save_path, 'probe2.txt'), 'w')
    x3_file = open(os.path.join(save_path, 'probe3.txt'), 'w')
    x4_file = open(os.path.join(save_path, 'probe4.txt'), 'w')
    x5_file = open(os.path.join(save_path, 'probe5.txt'), 'w')
    x6_file = open(os.path.join(save_path, 'probe6.txt'), 'w')
    x7_file = open(os.path.join(save_path, 'probe7.txt'), 'w')
    return x1_file, x2_file, x3_file, x4_file, x5_file, x6_file, x7_file

def save_probes(t, h_n0, dw_beach, x1_file, x2_file, x3_file, x4_file, x5_file, x6_file, x7_file, Ind_1, Ind_2, Ind_3, Ind_4, Ind_5, Ind_6, Ind_7):

    #------------------- wave elevation probe 1 : x1 = 15m -------------------#
    x1_file.write('%-10s %-10s\n'
                  %(str(t),str(h_n0.dat.data[Ind_1]+dw_beach.dat.data[Ind_1])))
    #------------------- wave elevation probe 2 : x2 = 17m -------------------#
    x2_file.write('%-10s %-10s\n'
                  %(str(t),str(h_n0.dat.data[Ind_2]+dw_beach.dat.data[Ind_2])))
    #------------------- wave elevation probe 3 : x3 = 19m -------------------#
    x3_file.write('%-10s %-10s\n'
                  %(str(t),str(h_n0.dat.data[Ind_3]+dw_beach.dat.data[Ind_3])))
    #------------------- wave elevation probe 4 : x4 = 20m -------------------#
    x4_file.write('%-10s %-10s\n'
                  %(str(t),str(h_n0.dat.data[Ind_4]+dw_beach.dat.data[Ind_4])))
    #------------------- wave elevation probe 5 : x5 = 21m -------------------#
    x5_file.write('%-10s %-10s\n'
                  %(str(t),str(h_n0.dat.data[Ind_5]+dw_beach.dat.data[Ind_5])))
    #------------------- wave elevation probe 6 : x6 = 22m -------------------#
    x6_file.write('%-10s %-10s\n'
                  %(str(t),str(h_n0.dat.data[Ind_6]+dw_beach.dat.data[Ind_6])))
    #------------------- wave elevation probe 7 : x7 = 23m -------------------#
    x7_file.write('%-10s %-10s\n'
                  %(str(t),str(h_n0.dat.data[Ind_7]+dw_beach.dat.data[Ind_7])))

""" **************************************
    *   Updated interpolation with time  *
    **************************************"""
def update_input_data(ind, t, t_data, wm_data, \
                      wm_vel_data, dWM_dt_expr, \
                      WM_expr):
    #  if t<t2, update the interpolated interval:
    if t_data[ind+1]<=t:           # if t>t2:
        ind +=1
        #--------- update t1 and t2 --------#
        dWM_dt_expr.t1 = t_data[ind]
        dWM_dt_expr.t2 = t_data[ind+1]
        WM_expr.t1 = t_data[ind]
        WM_expr.t2 = t_data[ind+1]
        #----- update the velocity data ----#
        dWM_dt_expr.dwm1 = wm_vel_data[ind]
        dWM_dt_expr.dwm2 = wm_vel_data[ind+1]
        #------ update the motion data -----#
        WM_expr.wm1 = wm_data[ind]
        WM_expr.wm2 = wm_data[ind+1]
    
    return ind, dWM_dt_expr, WM_expr

