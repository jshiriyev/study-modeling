"""This is going to be all about finite difference!"""

from scipy.sparse import csr_matrix as csr
from scipy.sparse import csc_matrix as csc
from scipy.sparse import hstack
from scipy.sparse import triu

from scipy.sparse.linalg import spsolve_triangular as sps
from scipy.sparse.linalg import splu

import matplotlib.pyplot as plt
import numpy as np

# from mesh import rectangle

class Poisson():

    """
    Finite Difference Solution of 2D Poisson's equation

    # \\del^2 U/\\del x^2+\\del^2 U/\\del y^2 = f(x,y)

    with boundary conditions:

    u = 0 at left and right edges,
    u = 10 at bottom edge,
    and at top edge U+\\del U/\\del y = 5

    """

    def __init__(self,rectangle):

        self.rectangle = rectangle

    def solve(self):

        pass

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