"""
Created on Fri Oct 23 14:10:14 2015

@author: mmtjs

Class for the beam.
"""

import firedrake as fd
import numpy as np

class Beam(object):
    """
    Describes beam domain, FE discretization and solution.    
    """
    def __init__(self, nx, ny, nz, Lx, Ly, Lz, rho, lam, mu, dt, t_end, **kwargs):
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.dt = dt
        self.t_end = t_end
        self.rho = rho
        self.lam = lam
        self.mu = mu
        self.t = np.arange(0., self.t_end + self.dt, self.dt)
        self.nonlin = kwargs['nonlin']
        self.g = kwargs['g']
        self.initialize()

    def initialize(self):
        self.bottom_id = 5
        self.interface_id = 1
        self.build_mesh()
        self.build_function_spaces()
        self.define_functions()
        self.initial_conditions()
        self.initialize_solvers()
        self.set_output_files()
        
    def build_mesh(self):
        self.mesh = fd.BoxMesh(self.nx, self.ny, self.nz, self.Lx, self.Ly, self.Lz)
#        self.mesh = fd.Mesh("beam.msh")
        
    def shift_mesh_origin(self, x, y):
        self.mesh.coordinates.dat.data[:,0] += x
        self.mesh.coordinates.dat.data[:,1] += y
           
    def build_function_spaces(self):
        self.V = fd.VectorFunctionSpace(self.mesh, "CG", 1)
        
    def define_functions(self):
        self.X = fd.Function(self.V, name="X") # displacement, not position
        self.U = fd.Function(self.V, name="U") # velocity
        self.trial = fd.TrialFunction(self.V)
#        self.phi_vect = fd.Function(self.V, name="Surface forcing")
        #Fgr = Function(V, name="Gravity")
        self.v = fd.TestFunction(self.V)
        
    def initial_conditions(self):
#  old:      ic1 = fd.project( fd.Expression([0.,0.,0.]), self.V )
#  old:      ic2 = fd.project( fd.Expression(["0.1*(1-cos(pi*x[2]/Lz/2.))",0.,0.], Lz=self.Lz), self.V )
        ic1 = fd.project( fd.Constant([0.,0.,0.]), self.V )
        x = fd.SpatialCoordinate(self.mesh)
        ic2 = fd.project( fd.as_vector([0.1*(1-fd.cos(fd.pi*x[2]/self.Lz/2.)),0.,0.]), self.V )
#        ic3 = fd.project( fd.Expression(["0.1*x[2]/Lz_B", 0., 0.], Lz_B=self.Lz), self.V )
        self.X.assign(ic2)
#        self.static_solver()
        self.U.assign(ic1)
#        self.read_raw()

    def surface_force(self):
        A = 1.
        d = 0.9*self.Lz
        l = 0.1*self.Lz
        T = fd.Function(self.V)
        T.interpolate( fd.conditional( self.X[2]>d,fd.as_vector([A*fd.sin((self.X[2]-d)/l * fd.pi/2.)**2,0.,0.]), fd.as_vector([0., 0., 0.]) ) )  # surface force / area
        return T   

# wrong it is the variable X not the coordinate x: x = fd.SpatialCoordinate(self.mesh)
# wrong def eval(self, value, X):

#  old: class MyExpression(fd.Expression):
#  old:     def eval(self, value, X):
#  old:         if X[2]>d:
#  old:                    value[:] = [A*fd.sin((X[2]-d)/l * fd.pi/2.)**2,0.,0.]
#  old:                else:
#  old:                    value[:] = [0.,0.,0.]
#  old:     def value_shape(self):
#  old:         return (3,)
#  old:  T.interpolate(MyExpression()) # surface force / areaT
#  old:  return T   
        
    def static_solver(self):
        def epsilon(u):
            return 0.5*(fd.nabla_grad(u) + fd.nabla_grad(u).T)
            #return sym(nabla_grad(u))
        
        def sigma(u):
            d = u.geometric_dimension()  # space dimension
            return self.lam*fd.nabla_div(u)*fd.Identity(d) + 2*self.mu*epsilon(u)
        
        # Define variational problem
        u = fd.TrialFunction(self.V)
        v = fd.TestFunction(self.V)
#        f = fd.Constant((0, 0, -self.g)) # body force / rho
        f = fd.Constant((0, 0, 0)) # body force / rho
        T = self.surface_force()

        a = fd.inner(sigma(u), epsilon(v))*fd.dx
        L = fd.dot(f, v)*fd.dx + fd.dot(T, v)*fd.ds(1)
        # Compute solution
#  old:        self.DBC = fd.DirichletBC( self.V, fd.Expression([0.,0.,0.]), self.bottom_id)
        self.DBC = fd.DirichletBC( self.V, fd.as_vector([0.,0.,0.]), self.bottom_id)
        fd.solve(a == L, self.X, bcs=self.DBC)
        
     
        
    def initialize_solvers(self):
        # Kinematics                # Right Cauchy-Green tensor
        if self.nonlin:
            d = self.X.geometric_dimension()
            I = fd.Identity(d)             # Identity tensor
            F = I + fd.grad(self.X)             # Deformation gradient
            C = F.T*F   
            E = (C-I)/2.               # Green-Lagrangian strain
#            E = 1./2.*( fd.grad(self.X).T + fd.grad(self.X) + fd.grad(self.X).T * fd.grad(self.X) ) # alternative equivalent definition
        else:
            E = 1./2.*( fd.grad(self.X).T + fd.grad(self.X) ) # linear strain
            
        self.W = (self.lam/2.)*(fd.tr(E))**2 + self.mu*fd.tr( E*E )
#        f = fd.Constant((0, 0, -self.g)) # body force / rho
#        T = self.surface_force()
        
        # Total potential energy
        Pi = self.W * fd.dx
        # Compute first variation of Pi (directional derivative about X in the direction of v)
        F_expr = fd.derivative(Pi, self.X, self.v) 
       
        self.DBC = fd.DirichletBC( self.V, fd.as_vector([0.,0.,0.]), self.bottom_id)
        
#        delX = fd.nabla_grad(self.X)
#        delv_B = fd.nabla_grad(self.v)
#        T_x_dv = self.lam * fd.div(self.X) * fd.div(self.v) \
#                + self.mu * ( fd.inner( delX, delv_B + fd.transpose(delv_B) ) )
        
        self.a_U = fd.dot( self.trial, self.v ) * fd.dx
#        self.L_U = ( fd.dot( self.U, self.v ) - self.dt/2./self.rho * T_x_dv ) * fd.dx
        self.L_U = fd.dot( self.U, self.v ) * fd.dx - self.dt/2./self.rho * F_expr #\
#                  + self.dt/2./self.rho*fd.dot(T,self.v)*fd.ds(1) # surface force at x==0 plane
                  # + self.dt/2.*fd.dot(f,self.v)*fd.dx # body force
        self.a_X = fd.dot( self.trial, self.v ) * fd.dx
#        self.L_interface = fd.dot(self.phi_vect, self.v) * fd.ds(self.interface_id)
        self.L_X = fd.dot( (self.X + self.dt * self.U), self.v ) * fd.dx #\
#                    - self.dt/self.rho * self.L_interface
        
        self.LVP_U = fd.LinearVariationalProblem(self.a_U, self.L_U, self.U, bcs = [self.DBC])
        self.LVS_U = fd.LinearVariationalSolver(self.LVP_U)
        self.LVP_X = fd.LinearVariationalProblem(self.a_X, self.L_X, self.X, bcs = [self.DBC])
        self.LVS_X = fd.LinearVariationalSolver(self.LVP_X)
    
    def set_output_files(self):
        self.outfile_X = fd.File("results/X.pvd")
        self.outfile_U = fd.File("results/U.pvd")
        
    def output_data(self):
        mesh_static = self.mesh.coordinates.vector().get_local()
        self.mesh.coordinates.vector().set_local( mesh_static + self.X.vector().get_local() )
        self.outfile_X.write( self.X )
        self.outfile_U.write( self.U )
        self.mesh.coordinates.vector().set_local( mesh_static )
        
    def write_raw(self):
        rawfile = fd.DumbCheckpoint("X_end", mode=fd.FILE_CREATE)
        rawfile.store(self.X)
        rawfile.close()
        rawfile = fd.DumbCheckpoint("U_end", mode=fd.FILE_CREATE)
        rawfile.store(self.U)
        rawfile.close()
        
    def read_raw(self):
        rawfile = fd.DumbCheckpoint("X_end", mode=fd.FILE_READ)
        rawfile.load(self.X)
        rawfile.close()
        rawfile = fd.DumbCheckpoint("U_end", mode=fd.FILE_READ)
        rawfile.load(self.U)
        rawfile.close()       
        
    def evolve_time(self):
        # Stoermer-Verlet scheme
        self.LVS_U.solve()
        self.LVS_X.solve()
        self.LVS_U.solve()
        
    def E_pot(self):
#        delX = fd.nabla_grad(self.X)
#        divX = fd.div(self.X)
#        return fd.assemble( 1./2. * ( self.lam*(divX**2) + self.mu*( fd.inner( delX, delX + fd.transpose(delX) ) ) ) * fd.dx )
        return fd.assemble( self.W * fd.dx )
        
    def E_kin(self):
        return fd.assemble( 1./2. * self.rho * fd.dot(self.U, self.U) * fd.dx )
