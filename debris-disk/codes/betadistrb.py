"""
Functions drawing radiation beta 
from different beta distributions

@author: Eve J. Lee
May 17th 2016
"""

import numpy as np
import numpy.random as nr
import scipy.integrate as sint
import scipy.interpolate as si
import consts
 
def Donhanyi(e, cosf, betapow=1.5, betamin=0.001, betamax=1, Ndust=100000, beta_bounded=False, a=50, Mstar=1, Tage=float('inf')):
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1./3.)
        f_betamax = (1. - e ** 2) * (a - half_period_term)/(a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf))
    else: 
        f_betamax = (1-e**2)/2./(1+e*cosf)
    if f_betamax > 1: f_betamax = betamax
    Pmax = (f_betamax**(1+betapow)-betamin**(1+betapow))/(betamax**(1+betapow)-betamin**(1+betapow))
    return (nr.uniform(0, Pmax, Ndust)*(betamax**(1+betapow)-betamin**(1+betapow))+betamin**(1+betapow))**(1./(1+betapow))

def f_Donhanyi(betapow=1.5, betamin=0.001, betamax=1):
    f = lambda P: (P*(betamax**(1+betapow)-betamin**(1+betapow))+betamin**(1+betapow))**(1./(1+betapow))
    return np.vectorize(f)

def OrbTimeCorr(e, cosf, betapow=1.5, betamin=0.001, betamax=1, Ndust=100000, stabfac = 0.998, beta_bounded=False, a=50, Mstar=1, Tage=float('inf')):
    dNdbeta = lambda beta: beta**betapow*(1-beta)**1.5*(1-e**2-2*beta*(1+e*cosf))**-1.5
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        f_betamax = stabfac * (1. - e ** 2) * (a - half_period_term) / (
                    a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf))
    else:
        f_betamax = (1 - e ** 2) / 2. / (1 + e * cosf) * stabfac
    if f_betamax > 1: f_betamax = betamax*stabfac
    norm = sint.quad(dNdbeta, betamin, f_betamax)[0]

    Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0]/norm
    Nbeta = np.vectorize(Nbeta)

    invNbeta = si.interp1d(Nbeta(np.linspace(betamin*0.8, f_betamax, 20)), np.linspace(betamin*0.8, f_betamax, 20))
    return invNbeta(nr.uniform(0, 1, Ndust))

def OrbTimeCorrCirc(betapow=1.5, betamin=0.001, Ndust=100000, stabfac = 0.998, beta_bounded=False, a=50, Mstar=1, Tage=float('inf')):
    dNdbeta = lambda beta: beta**betapow*((1-beta)/(1-2*beta))**1.5
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        f_betamax = (half_period_term - a) / (2 * half_period_term - a) * stabfac
    else:
        f_betamax = 0.5*stabfac

    if f_betamax > 0.5:
        f_betamax = 0.5 * stabfac
    norm = sint.quad(dNdbeta, betamin, f_betamax)[0]

    Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0]/norm
    Nbeta = np.vectorize(Nbeta)

    invNbeta = si.interp1d(Nbeta(np.linspace(betamin*0.8, f_betamax, 20)), np.linspace(betamin*0.8, f_betamax, 20))
    
    return invNbeta(nr.uniform(0, 1, Ndust))

