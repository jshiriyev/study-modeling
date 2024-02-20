"""
2500 Solved Problems in Fluid Mechanics and Hydraulics

by Jack B. Evett and Cheng Liu

Problem 7.2: What angle \alpha of jet is required to reach
the roof of the building with minimum jet velocity at the
nozzle? What is the value of initial velocity?

H = 28 m (height of the roof of the building)
L = 24 m (lateral distance to the building)

"""

from matplotlib import pyplot as plt

import numpy as np

def projectile_motion(vo,alpha):

    tmax = vo*np.sin(np.deg2rad(alpha))/9.81

    time = np.linspace(0,2*tmax,1000)

    x = vo*np.cos(np.deg2rad(alpha))*time
    y = vo*np.sin(np.deg2rad(alpha))*time-9.81*time**2/2

    return x,y

plt.plot(*projectile_motion(25.2,69.7),linestyle='--')

plt.vlines(24,0,28,color='red')

plt.ylim((0,35))

plt.xlabel("Lateral Distance")
plt.ylabel("Vertical Distance")

plt.show()
