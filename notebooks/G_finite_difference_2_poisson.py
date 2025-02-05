import matplotlib.pyplot as plt

import numpy as np

from scipy.sparse import csr_matrix as csr
from scipy.sparse import csc_matrix as csc
from scipy.sparse import hstack
from scipy.sparse import triu

from scipy.sparse.linalg import spsolve_triangular as sps
from scipy.sparse.linalg import splu

def poisson1D(self):
    
    """Start Input"""
    xN0 = 0
    xN1 = 7
    beta = 10
    bL = (1,0,200)
    bU = (0,1,50)
    """End Input"""

    grids = computational_singlephase()

    grids.cartesian((xN1,1,1),(70,1,1))

    grids.initialize()

    grids.transmissibility()

    grids.central(order=2)

    grids.implement_bc(bL,bU)

    grids.solve(rhs=beta)

    ##analytical = diffusivity()
    ##
    ##analytical.cartesian_poisson_1D((xN0,xN1),(bL,bU),beta)
    
    ##plt.scatter(grids.center[:,0],grids.unknown,c='r',label='Finite Difference')
    ##
    ####plt.plot(analytical.x,analytical.u,'k',label='Analytical')
    ##
    ##plt.xlabel('x-axis')
    ##plt.ylabel('pressure')
    ##
    ##plt.legend(bbox_to_anchor=(0., 1.02, 1., .102),
    ##           loc='lower left',ncol=2, mode="expand", borderaxespad=0.)
    ##
    ##plt.xlim([0,7])
    ####plt.ylim([0,None])
    ##
    ##plt.show()

def poisson2D(self):

    solver = computational_singlephase()

    num_x = 5
    num_y = 5

    solver.cartesian((1.,1.,1.),(num_x,num_y,1))

    solver.initialize()

    solver.transmissibility()

    solver.central()

    solver.implement_bc(b_xmin=(1,0,0),
                        b_xmax=(1,0,0),
                        b_ymin=(1,0,10),
                        b_ymax=(1,1,5))

    solver.solve(rhs=solver.center[:,0]**2)

    X = solver.center[:,0].reshape(num_x,num_y)
    Y = solver.center[:,1].reshape(num_x,num_y)

    P = solver.unknown.reshape(num_x,num_y)

    ##plt.contourf(X,Y,P,alpha=1,cmap=cm.PuBu)
    ##plt.colorbar()
    ##
    ##plt.title('Pressure Map',fontsize=14)
    ##plt.xlabel('x-axis',fontsize=14)
    ##plt.ylabel('y-axis',fontsize=14)
    ##
    ##plt.xlim([0,1])
    ##plt.ylim([0,1])
    ##
    ##plt.show()

class Poisson():

    """
    Finite Difference Solution of 2D Poisson's equation

    \\del^2 U/\\del x^2+\\del^2 U/\\del y^2 = f(x,y)

    with boundary conditions:

    u = 0 at left and right edges,
    u = 10 at bottom edge,
    and at top edge U+\\del U/\\del y = 5

    """

    def __init__(self,rectangle):

        self.rectangle = rectangle

    def solve(self):

        pass

if __name__ == "__main__"

    Nx = 5
    Ny = 5

    dx = 1./Nx
    dy = 1./Ny

    x = np.tile(np.linspace(-0.5+dx,0.5-dx,Nx-1),Ny)

    y = np.repeat(np.linspace(0.5,-0.5+dy,Ny),Nx-1)

    A1 = np.zeros(((Nx-1)*Ny,(Nx-1)*Ny))
    A2 = np.zeros(((Nx-1)*Ny,(Nx-1)*Ny))

    b = x**2

    ##xmin = np.array([0,4,8,12,16])
    xmin = np.arange(Ny)*(Nx-1)
    ##xmax = np.array([3,7,11,15,19])
    xmax = np.arange(Nx-2,(Nx-1)*Ny,Nx-1)

    ##ymin = np.array([16,17,18,19])
    ymin = np.arange((Nx-1)*(Ny-1),(Nx-1)*Ny)
    ##ymax = np.array([0,1,2,3])
    ymax = np.arange(0,Nx-1)

    for i in range((Nx-1)*Ny):
        
        A1[i,i] = -2/dx**2
        
        if not np.any(i==xmin):
            A1[i,i-1] = 1/dx**2
            
        if not np.any(i==xmax):
            A1[i,i+1] = 1/dx**2
            
        A2[i,i] = -2/dy**2
        
        if not np.any(i==ymin):
            A2[i,i+Nx-1] = 1/dy**2
        else:
            b[i] -= 10/dy**2

        if not np.any(i==ymax):
            A2[i,i-Nx+1] = 1/dy**2
        else:
            A2[i,i] -= 2/dy
            A2[i,i+Nx-1] = 2/dy**2
            b[i] -= 10/dy

    A = A1+A2

    # print(b)

    u = np.linalg.solve(A,b)

    print(u)

    # X = x.reshape((Ny,Nx-1))
    # Y = y.reshape((Ny,Nx-1))
    # U = u.reshape((Ny,Nx-1))

    # plt.contourf(X,Y,U,cmap='hot')
    # plt.xlim((-0.5,0.5))
    # plt.ylim((-0.5,0.5))
    # plt.colorbar()
    # plt.show()

    """
    LU decomposition

    """

    A = csc(A)

    u_LU = splu(A).solve(b)

    print(u_LU)


    """
    Gauss Elimination Method
    """

    # Forward elimination

    # a_ij = a_ij - a_ik/a_kk*a_kj
    # k: diagonals -- 0 to n-2
    # i: rows -- k+1 to n-1
    # j: cols -- k to n-1

    # print(augmented_upper[:,0:-1].toarray().shape)

    # u_GE = sps(augmented_upper[:,0:-1],augmented_upper[:,-1].toarray().flatten(),lower=False)

    # print(augmented_upper[:,0:-1].toarray())
    # print(augmented_upper[:,-1].toarray().flatten())

    # print(u_GE)

    """
    Gauss-Seidel Method is an iterative approach
    """