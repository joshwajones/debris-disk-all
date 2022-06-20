"""
Functions to draw particles/wires

@author: Eve J. Lee
Feb. 15th 2016
"""

import numpy as np
import prob_position_orbit as ppo
import pdb

def draw_wire(ax, omega, I, argu_peri = 0, e = 0, r0=1, format='b-', lw=2): #draw particles in orbits
    theta = np.linspace(-2 * np.pi, 2 * np.pi, 100)
    r = r0*(1-e**2.)/(1+e*np.cos(theta))
    x = r*(np.cos(omega)*np.cos(theta+argu_peri)-\
        np.sin(omega)*np.sin(theta+argu_peri)*np.cos(I))
    y = r*(np.sin(omega)*np.cos(theta+argu_peri)+\
        np.cos(omega)*np.sin(theta+argu_peri)*np.cos(I))
    z = r*np.sin(argu_peri+theta)*np.sin(I)
    ax.plot(x, y, z, format, lw=lw)

def draw_particles(ax, omega, I, argu_peri = 0, e = 0, r0=1, Npts=100, format='bo', \
                   msz=0.5, ylim=(-20,20), zlim=(-10,10)): #draw particles in orbits
    theta, ftheta = ppo.OutputPosition(e, Npts=Npts)
    r = r0*(1-e**2.)/(1+e*np.cos(theta))
    x = r*(np.cos(omega)*np.cos(theta+argu_peri)-\
        np.sin(omega)*np.sin(theta+argu_peri)*np.cos(I))
    y = r*(np.sin(omega)*np.cos(theta+argu_peri)+\
        np.cos(omega)*np.sin(theta+argu_peri)*np.cos(I))
    z = r*np.sin(argu_peri+theta)*np.sin(I)
    ax.plot(x, y, 0, format, markersize=msz, markeredgecolor='none')
