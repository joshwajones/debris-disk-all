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
import math
import pdb 

rand_number = 0.7991
 
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

def OrbTimeCorr(seed=-1, e=0, cosf=0, betapow=1.5, betamin=0.001, betamax=1, Ndust=100000, stabfac = 0.998, beta_bounded=False, a=50, Mstar=1, Tage=float('inf')):
    dNdbeta = lambda beta: beta**betapow*(1-beta)**1.5*(1-e**2-2*beta*(1+e*cosf))**-1.5
    f_betamax = (1 - e ** 2) / 2. / (1 + e * cosf) * stabfac
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        maxbeta_calc = stabfac * (1. - e ** 2) * (a - half_period_term) / (
                    a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf))
        maxbeta_calc = max(0, maxbeta_calc)
        f_betamax = min(f_betamax, maxbeta_calc)
    if f_betamax > 1: f_betamax = betamax*stabfac
    
    norm = sint.quad(dNdbeta, betamin, f_betamax)[0]

    Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0]/norm
    Nbeta = np.vectorize(Nbeta)

    invNbeta = si.interp1d(Nbeta(np.linspace(betamin*0.8, f_betamax, 20)), np.linspace(betamin*0.8, f_betamax, 20))
    # 0.7991
    if seed < 0:
        seed = nr.uniform(0, 1, Ndust)
    return  invNbeta(seed), f_betamax
    #return invNbeta(nr.uniform(0, 1, Ndust)), f_betamax

def OrbTimeCorr_legacy(e, cosf, betapow=1.5, betamin=0.001, betamax=1, Ndust=100000, stabfac = 0.998, beta_bounded=False, a=50, Mstar=1, Tage=float('inf')):
    dNdbeta = lambda beta: beta**betapow*(1-beta)**1.5*(1-e**2-2*beta*(1+e*cosf))**-1.5
    f_betamax = (1 - e ** 2) / 2. / (1 + e * cosf) * stabfac
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        maxbeta_calc = stabfac * (1. - e ** 2) * (a - half_period_term) / (
                a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf))
        maxbeta_calc = max(0, maxbeta_calc)
        f_betamax = min(f_betamax, maxbeta_calc)
    if f_betamax > 1: f_betamax = betamax * stabfac
    norm = sint.quad(dNdbeta, betamin, f_betamax)[0]

    Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0]/norm
    Nbeta = np.vectorize(Nbeta)

    invNbeta = si.interp1d(Nbeta(np.linspace(betamin*0.8, f_betamax, 20)), np.linspace(betamin*0.8, f_betamax, 20))
    return invNbeta(nr.uniform(0, 1, Ndust))



def OrbTimeCorr_opt(e, cosf, betapow=1.5, betamin=0.001, betamax=0.5, Ndust=100000, stabfac = 0.998, beta_bounded=False, a=50, Mstar=1, Tage=float('inf')):
    a *= consts.au2cm
    mu = consts.G * Mstar
    def get_betamax():
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        f_betamax = stabfac * (1. - e ** 2) * (a - half_period_term) / (
                a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf))
        print(f_betamax)
    Tage = 1
    for _ in range(100):
        get_betamax()
        Tage *= 2
    print((1 - e ** 2) / 2. / (1 + e * cosf) * stabfac)
    pdb.set_trace()
    #approximate e~0
    eprox = 0.
    dNdbeta = lambda beta: beta**betapow*(1-beta)**1.5*(1-e**2-2*beta*(1+e*cosf))**-1.5
    dNdbeta2 = lambda beta: beta**betapow*(1-beta)**1.5*(1-eprox**2-2*beta*(1+eprox*cosf))**-1.5
    f_betamax = (1 - e ** 2) / 2. / (1 + e * cosf) * stabfac
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        f_betamax = min(f_betamax, stabfac * (1. - e ** 2) * (a - half_period_term) / (
                a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf)))

    if f_betamax > 1: f_betamax = betamax*stabfac
    temp = e
    e = 0
    f_betaprox = (1 - e ** 2) / 2. / (1 + e * cosf) * stabfac
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        f_betaprox = min(f_betaprox, stabfac * (1. - e ** 2) * (a - half_period_term) / (
                a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf)))

    #f_betaprox = 0.5 * stabfac
    print(f_betaprox)
    e = temp



    norm = sint.quad(dNdbeta, betamin, f_betamax)[0]
    normprox = sint.quad(dNdbeta2, betamin, f_betaprox)[0]

    Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0]/norm
    Nbeta = np.vectorize(Nbeta)

    Nbeta2 = lambda beta: sint.quad(dNdbeta2, betamin, beta)[0]
    Nbeta2 = np.vectorize(Nbeta2)

    #inverse = si.interp1d(Nbeta2(np.linspace(betamin*0.8, f_betamax, 20)), np.linspace(betamin*0.8, f_betamax, 20))
    invNbeta = si.interp1d(Nbeta(np.linspace(betamin*0.8, f_betamax, 20)), np.linspace(betamin*0.8, f_betamax, 20))
    #inverse4 = si.interp1d(Nbeta2(np.linspace(betamin*0.8, f_betamax, 20)), np.linspace(betamin*0.8, f_betamax, 20))
    inverse3 = si.interp1d(Nbeta2(np.linspace(betamin*0.8, f_betaprox, 20)), np.linspace(betamin*0.8, f_betaprox, 20), bounds_error=False, fill_value="extrapolate")
    vals = nr.uniform(0, 1, Ndust)
    # print(norm)
    # print(betamin*0.8, f_betamax)
    #print(invNbeta(vals))

    #print(inverse(vals * norm))

    trouble = inverse3(vals * normprox)
    # if math.isnan(trouble):
    #     print("val:  ", vals)
    #     pdb.set_trace()
    # print(trouble)
    # tr = inverse4(vals * norm)
    # if math.isnan(tr):
    #     pdb.set_trace()
    # print(tr)
    # print()
    # print()
    # print()
    print("yo")
    print(f_betaprox)
    print(f_betamax)
    print(invNbeta(vals), trouble)
    return invNbeta(vals), trouble

def get_inverse_CDF(e=0, cosf=1, betapow=1.5, betamin=0.001, betamax=1, Ndust=100000, stabfac = 0.998, beta_bounded=False, a=50, Mstar=1, Tage=float('inf'), precision=100):
    dNdbeta = lambda beta: beta ** betapow * (1 - beta) ** 1.5 * (1 - e ** 2 - 2 * beta * (1 + e * cosf)) ** -1.5
    f_betamax = (1 - e ** 2) / 2. / (1 + e * cosf) * stabfac
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        maxbeta_calc = stabfac * (1. - e ** 2) * (a - half_period_term) / (
                a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf))
        maxbeta_calc = max(maxbeta_calc, 0) #Handle case where < 0? happens for small-medium ages
        f_betamax = min(f_betamax, maxbeta_calc)
    if f_betamax > betamax:
        f_betamax = stabfac * betamax

    norm = sint.quad(dNdbeta, betamin, f_betamax)[0]
    Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0]/norm
    Nbeta = np.vectorize(Nbeta)
    invNbeta = si.interp1d(Nbeta(np.linspace(betamin * 0.8, f_betamax, precision)), np.linspace(betamin * 0.8, f_betamax, precision))
    
    return invNbeta, f_betamax


def get_approx_interp(e, cosf, betapow=1.5, betamin=0.001, betamax=1, Ndust=100000, stabfac = 0.998, beta_bounded=False, a=50, Mstar=1, Tage=float('inf'), precision=100):


    e = 0 # approximate e~0
    dNdbeta = lambda beta: beta ** betapow * (1 - beta) ** 1.5 * (1 - e ** 2 - 2 * beta * (1 + e * cosf)) ** -1.5
    f_betamax = (1 - e ** 2) / 2. / (1 + e * cosf) * stabfac
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        maxbeta_calc = stabfac * (1. - e ** 2) * (a - half_period_term) / (
                a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf))
        maxbeta_calc = max(maxbeta_calc, 0)
        f_betamax = min(f_betamax, maxbeta_calc)
    if f_betamax > betamax: 
        f_betamax = stabfac * betamax
    norm = sint.quad(dNdbeta, betamin, f_betamax)[0]
    Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0]/norm
    Nbeta = np.vectorize(Nbeta)
    invNbeta_unscaled = si.interp1d(Nbeta(np.linspace(betamin * 0.8, f_betamax, precision)), np.linspace(betamin * 0.8, f_betamax, precision))
    return invNbeta_unscaled, f_betamax


def OrbTimeCorr_Vectorized(inverseCDF, Nlaunch):
    return inverseCDF(nr.uniform(0, 1, Nlaunch))


def OrbTimeCorr_opt_approx(seed, inverseCDF, f_betamax_approximate, e, cosf, betapow=1.5, betamin=0.001, betamax=1,
                           Ndust=100000, stabfac=0.998,
                           beta_bounded=False, a=50, Mstar=1, Tage=float('inf'), idx=0):
    f_betamax = (1 - e ** 2) / 2. / (1 + e * cosf) * stabfac
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        f_betamax = stabfac * (1. - e ** 2) * (a - half_period_term) / (
                a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf))
    if f_betamax > 1: f_betamax = betamax * stabfac

    # beta = inverseCDF(nr.uniform(0, 1, Ndust))
    beta = inverseCDF(seed)
    # exceeds_maxbeta = np.where(beta > f_betamax)
    # for i in range(len(exceeds_maxbeta)):
    #     print(beta, f_betamax)
    #     print((beta > f_betamax)[0])
    #     print(idx)
    #     pdb.set_trace()
    if (beta > f_betamax)[0]:
        # print(beta, f_betamax)
        # print((beta > f_betamax)[0])
        beta *= f_betamax / f_betamax_approximate * stabfac
        # print(idx)
        # pdb.set_trace()
    return beta

def OrbTimeCorr_opt_approx(seed, inverseCDF, f_betamax_approximate, e, cosf, betapow=1.5, betamin=0.001, betamax=1, Ndust=100000, stabfac=0.998,
                          beta_bounded=False, a=50, Mstar=1, Tage=float('inf'), idx = 0):
    
    f_betamax = (1 - e ** 2) / 2. / (1 + e * cosf) * stabfac
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        f_betamax = stabfac * (1. - e ** 2) * (a - half_period_term) / (
                a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf))
    if f_betamax > 1: f_betamax = betamax * stabfac
    
    #beta = inverseCDF(nr.uniform(0, 1, Ndust))
    beta = inverseCDF(seed)
    # exceeds_maxbeta = np.where(beta > f_betamax)
    # for i in range(len(exceeds_maxbeta)):
    #     print(beta, f_betamax)
    #     print((beta > f_betamax)[0])
    #     print(idx)
    #     pdb.set_trace()
    if (beta > f_betamax)[0]:
        # print(beta, f_betamax)
        # print((beta > f_betamax)[0])
        beta *= f_betamax / f_betamax_approximate * stabfac
        #print(idx)
        # pdb.set_trace()
    return beta

def OrbTimeCorr_optimized(inverseCDF, e, cosf, betapow=1.5, betamin=0.001, betamax=1, Ndust=100000, stabfac = 0.998, beta_bounded=False, a=50, Mstar=1, Tage=float('inf')):
    #approximate e~0
    eprox = 0.
    dNdbeta2 = lambda beta: beta**betapow*(1-beta)**1.5*(1-eprox**2-2*beta*(1+eprox*cosf))**-1.5
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        f_betamax = stabfac * (1. - e ** 2) * (a - half_period_term) / (
                    a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf))
    else:
        f_betamax = (1 - e ** 2) / 2. / (1 + e * cosf) * stabfac
    if f_betamax > 1: f_betamax = betamax*stabfac
    #f_betaprox = 0.5 * stabfac
    e = 0
    if beta_bounded:
        a *= consts.au2cm
        mu = consts.G * Mstar
        half_period_term = pow((mu * Tage ** 2 / (np.pi ** 2)), 1. / 3.)
        f_betaprox = stabfac * (1. - e ** 2) * (a - half_period_term) / (
                    a * (1. - e ** 2) - 2. * half_period_term * (1. + e * cosf))
    else:
        f_betaprox = (1 - e ** 2) / 2. / (1 + e * cosf) * stabfac
    

    norm = sint.quad(dNdbeta, betamin, f_betamax)[0]
    normprox = sint.quad(dNdbeta2, betamin, f_betaprox)[0]

    Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0]/norm
    Nbeta = np.vectorize(Nbeta)

    Nbeta2 = lambda beta: sint.quad(dNdbeta2, betamin, beta)[0]
    Nbeta2 = np.vectorize(Nbeta2)

    #inverse = si.interp1d(Nbeta2(np.linspace(betamin*0.8, f_betamax, 20)), np.linspace(betamin*0.8, f_betamax, 20))
    invNbeta = si.interp1d(Nbeta(np.linspace(betamin*0.8, f_betamax, 20)), np.linspace(betamin*0.8, f_betamax, 20))
    #inverse4 = si.interp1d(Nbeta2(np.linspace(betamin*0.8, f_betamax, 20)), np.linspace(betamin*0.8, f_betamax, 20))
    inverse3 = si.interp1d(Nbeta2(np.linspace(betamin*0.8, f_betaprox, 20)), np.linspace(betamin*0.8, f_betaprox, 20), bounds_error=False, fill_value="extrapolate")
    vals = nr.uniform(0, 1, Ndust)
    # print(norm)
    # print(betamin*0.8, f_betamax)
    #print(invNbeta(vals))

    #print(inverse(vals * norm))

    trouble = inverse3(vals * normprox)
    # if math.isnan(trouble):
    #     print("val:  ", vals)
    #     pdb.set_trace()
    # print(trouble)
    # tr = inverse4(vals * norm)
    # if math.isnan(tr):
    #     pdb.set_trace()
    # print(tr)
    # print()
    # print()
    # print()
    return invNbeta(vals), trouble

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

