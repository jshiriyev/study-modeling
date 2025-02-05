import matplotlib.pyplot as plt

import numpy as np

from scipy.sparse import csr_matrix as csr
from scipy.sparse import csc_matrix as csc
from scipy.sparse import hstack
from scipy.sparse import triu

from scipy.sparse.linalg import spsolve_triangular as sps
from scipy.sparse.linalg import splu

def laplace1D(self):
    
    bnd1 = (1,-1,50)
    bnd2 = (1,0,100)

    grids = computational_singlephase()

    grids.cartesian((7,1,1),(7,1,1))

    grids.initialize()

    grids.transmissibility()

    grids.central(order=2)

    grids.implement_bc(b_xmin=bnd1,b_xmax=bnd2)

    ##analytical = core_singlephase()
    ##analytical.cartesian_poisson_1D((0,7),(bnd1,bnd2),0)

    grids.solve()

def laplace2D(self):

    grids = computational_singlephase()

    num_x = 50
    num_y = 50

    grids.cartesian((20,10,10),(num_x,num_y,1))

    grids.initialize()

    grids.transmissibility()

    grids.central(order=2)

    grids.implement_bc(b_xmin=(1,0,0),
                       b_xmax=(0,1,0),
                       b_ymin=(1,0,0),
                       b_ymax=(1,0,100))

    grids.solve()

    ##X = grids.center[:,0].reshape(num_x,num_y)
    ##Y = grids.center[:,1].reshape(num_x,num_y)
    ##
    ##Z = grids.unknown.reshape(num_x,num_y)
    ##
    ##plt.contourf(X,Y,Z,alpha=1,cmap=cm.PuBu)
    ##plt.colorbar()
    ##
    ##plt.title('Pressure Map',fontsize=14)
    ##plt.xlabel('x-axis',fontsize=14)
    ##plt.ylabel('y-axis',fontsize=14)
    ##
    ##plt.xlim([0,20])
    ##plt.ylim([0,10])
    ##
    ##plt.show()

def laplace3D(self):

    solver = computational_singlephase()

    num_x = 25
    num_y = 25
    num_z = 10

    solver.cartesian((num_x*1.,num_y*1.,num_z*1.),(num_x,num_y,num_z))

    solver.initialize()

    solver.transmissibility()

    solver.central()

    solver.implement_bc(b_xmin=(1,0,0),
                        b_xmax=(0,1,0),
                        b_ymin=(1,0,0),
                        b_ymax=(1,0,100),
                        b_zmin=(0,1,0),
                        b_zmax=(0,1,0))

    solver.solve()

    ##X = solver.center[:,0].reshape(num_x,num_y,num_z)
    ##Y = solver.center[:,1].reshape(num_x,num_y,num_z)
    ##Z = solver.center[:,2].reshape(num_x,num_y,num_z)
    ##
    ####P = solver.unknown.reshape(num_x,num_y,num_z)
    ##
    ##fig = plt.figure()
    ##
    ##ax = fig.add_subplot(111, projection='3d')
    ##
    ##scatter = ax.scatter(X,Y,Z,alpha=0.5,c=solver.unknown)
    ##
    ##fig.colorbar(scatter,shrink=0.5,aspect=5)
    ##
    ##ax.set_xlim3d([0,solver.length_x])
    ##ax.set_ylim3d([0,solver.length_y])
    ##ax.set_zlim3d([0,solver.length_z])
    ##
    ##ax.set_xlabel('x-axis')
    ##ax.set_ylabel('y-axis')
    ##ax.set_zlabel('z-axis')
    ##
    ##plt.show()

class Laplace():

    def set_externalBC(self,b_xmin=(0,1,0),b_xmax=(0,1,0),b_ymin=(0,1,0),b_ymax=(0,1,0),b_zmin=(0,1,0),b_zmax=(0,1,0)):

        """
        b_xmin,b_xmax,b_ymin,b_ymax,b_zmin and b_zmax have three entries:
        - dirichlet boundary condition coefficient,
        - neumann boundary condition coefficient,
        - function value of boundary condition
        Default is no flow boundary conditions at the exterior boundaries
        """

        sizes = self.PorRock.grid_sizes

        indices = self.PorRock.grid_indices

        self.b_correction = np.zeros(self.PorRock.grid_numtot)

        questbound = lambda x: True if x>1 else False

        if questbound(self.PorRock.grid_num[0]):

            xmin = indices[:,0]==indices[:,1]
            xmax = indices[:,0]==indices[:,2]

            id_xmin = indices[xmin,0]
            id_xmax = indices[xmax,0]

            dx_xmin = sizes[id_xmin,0]
            dx_xmax = sizes[id_xmax,0]

            tx_xmin = self.transmissibility[id_xmin,0]
            tx_xmax = self.transmissibility[id_xmax,1]

            bc_xmin = (2*tx_xmin*dx_xmin)/(b_xmin[0]*dx_xmin-2*b_xmin[1])
            bc_xmax = (2*tx_xmax*dx_xmax)/(b_xmax[0]*dx_xmax+2*b_xmax[1])
            
            self.Gmatrix[id_xmin,id_xmin] -= bc_xmin*b_xmin[0]
            self.Gmatrix[id_xmax,id_xmax] -= bc_xmax*b_xmax[0]
            
            self.b_correction[id_xmin] -= bc_xmin*b_xmin[2]
            self.b_correction[id_xmax] -= bc_xmax*b_xmax[2]

        if questbound(self.PorRock.grid_num[1]):

            ymin = indices[:,0]==indices[:,3]
            ymax = indices[:,0]==indices[:,4]

            id_ymin = indices[ymin,0]
            id_ymax = indices[ymax,0]

            dy_ymin = sizes[id_ymin,1]
            dy_ymax = sizes[id_ymax,1]

            ty_ymin = self.transmissibility[id_ymin,2]
            ty_ymax = self.transmissibility[id_ymax,3]

            bc_ymin = (2*ty_ymin*dy_ymin)/(b_ymin[0]*dy_ymin-2*b_ymin[1])
            bc_ymax = (2*ty_ymax*dy_ymax)/(b_ymax[0]*dy_ymax+2*b_ymax[1])
            
            self.Gmatrix[id_ymin,id_ymin] -= bc_ymin*b_ymin[0]
            self.Gmatrix[id_ymax,id_ymax] -= bc_ymax*b_ymax[0]
            
            self.b_correction[id_ymin] -= bc_ymin*b_ymin[2]
            self.b_correction[id_ymax] -= bc_ymax*b_ymax[2]

        try:
            grid_num_x = self.PorRock.grid_num[2]
        except IndexError:
            grid_num_x = 1
            
        if questbound(grid_num_x):

            zmin = indices[:,0]==indices[:,5]
            zmax = indices[:,0]==indices[:,6]

            id_zmin = indices[zmin,0]
            id_zmax = indices[zmax,0]

            dz_zmin = sizes[id_zmin,2]
            dz_zmax = sizes[id_zmax,2]

            tz_zmin = self.transmissibility[id_zmin,4]
            tz_zmax = self.transmissibility[id_zmax,5]

            bc_zmin = (2*tz_zmin*dz_zmin)/(b_zmin[0]*dz_zmin-2*b_zmin[1])
            bc_zmax = (2*tz_zmax*dz_zmax)/(b_zmax[0]*dz_zmax+2*b_zmax[1])

            self.Gmatrix[id_zmin,id_zmin] -= bc_zmin*b_zmin[0]
            self.Gmatrix[id_zmax,id_zmax] -= bc_zmax*b_zmax[0]
            
            self.b_correction[id_zmin] -= bc_zmin*b_zmin[2]
            self.b_correction[id_zmax] -= bc_zmax*b_zmax[2]