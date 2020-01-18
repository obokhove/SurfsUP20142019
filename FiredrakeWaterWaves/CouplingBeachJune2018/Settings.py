
from firedrake import *

"""
    *********************************************
    *                 Test case                 *
    *********************************************"""
def input():
    #________________ Kind of data ________________#
    input_data = "measurements"   # from experiments
    #input_data = "created"      # set the wavemaker
    #_________________ Test case __________________#
    if input_data == "created":
        test_case = "test"    # choose a folder name
    elif input_data == "measurements":
        test_case = "111"        # choose from table
    return input_data, test_case


"""
    *********************************************************************
    *                         Numerical domain                          *
    *********************************************************************"""
def domain(input_data):
    if input_data == "created":
        #_______________________________ Beach _______________________________#
        H0 = 1.0                                  # Depth at rest (flat bottom)
        xb = 3.0                                           # Start of the beach
        sb = 0.2                                           # Slope of the beach
        H_expr = Expression("H0-0.5*(1+copysign(1.0,x[0]-xb))*slope*(x[0]-xb)",
                            H0=H0,xb=xb, slope=sb)
        
        #_______________________________ Basin _______________________________#
        Hend = 0.2                              # Depth at the end of the beach
        Lx = xb +(H0-Hend)/sb                                     # Length in x
        Lw = 1.0                                       # End of the x-transform
        n_res = 2
        res_x = 5.0*10**(-n_res)                                 # x-resolution
        n_z = 8                                        # Order of the expansion

    #----------------------------------------------------------------------#
    #                   WHEN COMPARED TO TUD EXPERIMENTS:                  #
    #----------------------------------------------------------------------#
    elif input_data == "measurements":
        #________________________________ Beach _______________________________#
        H0 = 1.005                                 # Depth at rest (flat bottom)
        xb = 20.055                                        # Start of the beach
        sb = 0.1                                           # Slope of the beach
        H_expr = Expression("H0-0.5*(1+copysign(1.0,x[0]-xb))*slope*(x[0]-xb)",
                            H0=H0,xb=xb, slope=sb)

        #________________________________ Basin _______________________________#
        Hend = 0.2                              # Depth at the end of the beach
        Lx = xb +(H0-Hend)/sb                                     # Length in x
        Lw = 1.0                                       # End of the x-transform
        n_res = 2
        res_x = 5.0*10**(-n_res)                                 # x-resolution
        n_z = 8                                        # Order of the expansion

    return H0, xb, sb, H_expr, Hend, Lx, Lw, res_x, n_z



"""
    **************************************************************************
    *                                Wavemaker                               *
    **************************************************************************"""
def wavemaker(H0, Lw, t):
    #_____________________________ Characteristics _____________________________#
    g = 9.81                                             # Gravitational constant
    lamb = 2.0                                                       # Wavelength
    k = 2*pi/lamb                                                   # Wave number
    w = sqrt(g*k*tanh(k*H0))                                     # Wave frequency
    Tw = 2*pi/w                                                     # Wave period
    gamma = 0.02                                                 # Wave amplitude
    t_stop = 40.0*Tw                                  # When to stop the wavemaker
    
    #________________________________ Expression _______________________________#
    WM_expr = \
        Expression("-0.5*(1+copysign(1.0,Lw-x[0]))*A*cos(w*t)",A=gamma, Lw = Lw, w=w, t=t)

    #_____________________________ Time derivative _____________________________#
    dWM_dt_expr = \
        Expression("0.5*(1+copysign(1.0,Lw-x[0]))*A*w*sin(w*t)",A=gamma, Lw=Lw, w=w, t=t)

    return g, lamb, k, w, Tw, gamma, t_stop, WM_expr, dWM_dt_expr



"""
    ***********************************
    *               Time              *
    ***********************************"""
def set_time():
    T0 = 0.0                # Initial time
    Tend = 50.0               # Final time
    t = T0             # Temporal variable
    dt = 0.001                 # Time step
    dt_save = 0.02      # saving time step
    return T0, t, dt, Tend, dt_save

