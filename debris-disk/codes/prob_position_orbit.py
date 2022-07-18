"""
Output position 
of particles in a given orbit

@author: Eve J. Lee
Feb. 16th 2016
"""

import numpy as np
import numpy.random as nr
import scipy.interpolate as sint
import matplotlib.pyplot as pl
import pdb

def GetMeanToTrueAnomalyFunc(e, Npts=100000):
    E = np.linspace(-np.pi, np.pi, Npts)
    M = E - e * np.sin(E)
    finvM = sint.interp1d(M, E)  # put this elsewhere for optimization
    return finvM
    

def OutputPosition(e, Npts=100, verbose=False):
    
    E = np.linspace(-np.pi, np.pi, Npts)
    M = E - e*np.sin(E)
    finvM = sint.interp1d(M, E) #put this elsewhere for optimization

    # calculates the tangent of interpolated true anomalies
    # choosing the mean anomaly uniformly
    tanf = np.sqrt((1+e)/(1-e))*\
           np.tan(finvM(nr.uniform(-np.pi, np.pi, Npts)/2.))

    theta = 2*np.arctan(tanf)
    if verbose:
        print(theta)
   
    theta = np.where(theta < 0, 2*np.pi+theta, theta)
    if verbose:
        print(theta)
 
    hist = np.histogram(theta, bins=20, normed=True)
    fhist = sint.interp1d(hist[1][1:], hist[0],bounds_error=False, fill_value=np.mean([hist[0][0], hist[0][-1]]))
    return theta, fhist


