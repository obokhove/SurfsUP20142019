
from firedrake import *

"""
    *********************************************
    *                 Test case                 *
    *********************************************"""
def test_case():
    #________________ Kind of data ________________#
    #input_data = "measurements"  # from experiments
    input_data = "created"       # set the wavemaker
    #______________ Temporal scheme _______________#
    scheme = "SE"
    #  "SE": Symplectic-Euler ; "SV": Stormer-Verlet
    #__________________ Dimension _________________#
    dim = "2D"
    #"2D": R(t) and b(x); "3D": R(y,t) and/or b(x,y)
    # if input = measurements, the dim must be 2D.
    #______ Path and name of the saved files ______#
    save_path = 'data/'+scheme+'/'+dim+'/'
    return input_data, scheme, dim, save_path


"""
    *********************************************************************
    *                         Numerical domain                          *
    *********************************************************************"""
def domain():
    #______________________ Beach ______________________#
    H0 = 1.0                                  # Depth at rest (flat bottom)
    xb = 3.0                                           # Start of the beach
    sb = 0.2                                           # Slope of the beach
    H_expr = Expression("H0-0.5*(1+copysign(1.0,x[0]-xb))*slope*(x[0]-xb)",
                        H0=H0,xb=xb, slope=sb)
    
    #______________________ Basin ______________________#
    Hend = 0.5                              # Depth at the end of the beach
    Lx = xb +(H0-Hend)/sb                                     # Length in x
    Ly = 1.0                                                  # Length in y
    Lw = 1.0                                       # End of the x-transform
    res_x = 0.05                                             # x-resolution
    res_y = 0.2                                             # y-resolution
    n_z = 8                                        # Order of the expansion
    return H0, xb, sb, H_expr, Hend, Lx, Ly, Lw, res_x, res_y, n_z

"""
    **************************************************************************
    *                                Wavemaker                               *
    **************************************************************************"""
def wavemaker(dim, H0, Ly, Lw, t):
    #_____________________________ Characteristics _____________________________#
    g = 9.81                                             # Gravitational constant
    lamb = 2.0                                                       # Wavelength
    k = 2*pi/lamb                                                   # Wave number
    w = sqrt(g*k*tanh(k*H0))                                     # Wave frequency
    Tw = 2*pi/w                                                     # Wave period
    gamma = 0.02                                                 # Wave amplitude
    t_stop = 40.0*Tw                                  # When to stop the wavemaker
    
    #________________________________ Expression _______________________________#
    if dim == "2D":
        WM_expr = \
        Expression("-0.5*(1+copysign(1.0,Lw-x[0]))*A*cos(w*t)",A=gamma, Lw = Lw, w=w, t=t)
    elif dim == "3D":
        WM_expr = \
        Expression("-0.5*(1+copysign(1.0,Lw-x[0]))*A*(x[1]-0.5*Ly)/(0.5*Ly)*cos(w*t)",A=gamma, Ly=Ly, Lw = Lw, w=w, t=t)
               
    #_____________________________ Time derivative _____________________________#
    if dim == "2D":
        dWM_dt_expr = \
        Expression("0.5*(1+copysign(1.0,Lw-x[0]))*A*w*sin(w*t)",A=gamma, Lw=Lw, w=w, t=t)
    elif dim == "3D":
        dWM_dt_expr = \
        Expression("0.5*(1+copysign(1.0,Lw-x[0]))*A*w*(x[1]-0.5*Ly)/(0.5*Ly)*sin(w*t)",A=gamma, Ly=Ly, Lw=Lw, w=w, t=t)
    #______________________________ y-derivative _______________________________#
    if dim == "2D":
        dWM_dy_expr = Expression("0.0")
    elif dim == "3D":
        dWM_dy_expr = Expression("-0.5*(1+copysign(1.0,Lw-x[0]))*A*cos(w*t)/(0.5*Ly)",A=gamma, Ly=Ly, Lw=Lw, w=w, t=t)

    return g, lamb, k, w, Tw, gamma, t_stop, WM_expr, dWM_dt_expr, dWM_dy_expr



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

