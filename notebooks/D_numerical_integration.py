import types

import numpy as np

import matplotlib.pyplot as plt

from scipy.integrate import quad
from scipy.integrate import trapezoid
from scipy.integrate import simps

def midpoint(f,x=None,bounds=(0,1),Nsample=100):

    if x is None:
        x = np.linspace(*bounds,Nsample)

    xmids = (xbnds[1:]+xbnds[:-1])/2

    if isinstance(f,types.FunctionType):
        ymids = f(xmids)
    else:
        ymids = (f[1:]+f[:-1])/2

    xdiff = (x[1:]-x[:-1])

    Inum = np.sum(ymids*xdiff)

    return Inum

def simpson(f,x=None,bounds=(0,1),Nsample=100):

    if x is None:
        xbnds = np.linspace(*bounds,Nsample)
    else:
        xbnds = x

    xmids = (xbnds[1:]+xbnds[:-1])/2

    if isinstance(f,types.FunctionType):
        ybnds = f(xbnds)
        ymids = f(xmids)
    else:
        ybnds = f
        ymids = (f[1:]+f[:-1])/2

    xdiff = (xbnds[1:]-xbnds[:-1])

    I1 = ybnds[:-1]*xdiff/6
    I2 = ymids*4*xdiff/6
    I3 = ybnds[1:]*xdiff/6

    Inum = np.sum(I1)+np.sum(I2)+np.sum(I3)

    return Inum

def trapezoid(f,x=None,bounds=(0,1),Nsample=100):

    if x is None:
        x = np.linspace(*bounds,Nsample)

    if isinstance(f,types.FunctionType):
        y = f(x)
    else:
        y = f

    xdiff = (x[1:]-x[:-1])
    ymids = (y[1:]+y[:-1])/2

    Inum = np.sum(ymids*xdiff)

    return Inum

def integrand(x):
    return 1/(1+x**2)

x = np.linspace(1,2,100)

y = integrand(x)

Iq = quad(integrand,1,2)
It = trapezoid(y,x)
Is = simps(y,x)

Itm = mytrapezoid(integrand,bounds=(1,2))
Ism = mysimpson(integrand,bounds=(1,2))

Iactual = np.arctan(2)-np.arctan(1)

## Exercise 2

##x = np.array([0,1.1,2.5,3.1,4.3,5.6,5.9,7.2,7.8,7.9,8.3,8.9,9.3,9.6,10.])
##y = np.array([0.5,0.691,0.845,0.880,0.896,0.846,0.826,0.714,0.652,0.642,0.598,0.533,0.489,0.457,0.414])
##
##It = trapezoid(y,x)
##Is = simps(y,x)

print("Scipy Quad gives {}".format(Iq))
print("SciPy Trapezoid gives {} and myTrapezoid gives {}".format(It,Itm))
print("SciPy Simpson gives {} and mySimpson gives {}".format(Is,Ism))
print("True result is {}".format(Iactual))
