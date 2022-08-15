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

def OrbTimeCorr_MidOptimized_Background_2(e_set, num_f=10, n_beta_grid=50, betapow=1.5, betamin=0.001, stabfac=0.997):

    #return inverse CDF for e, cosf combo
    #hashmap: e_val to list of cosf values

    e_vals = e_set
    cosf = np.linspace(-1, 1, num_f)
    f_vals = np.arccos(cosf)
    for i in range(len(f_vals)):
        if nr.uniform(0, 1) > 0.5:
            f_vals[i] = 2 * np.pi - f_vals[i]
    cosf = np.cos(f_vals)
    how_close = 1.e-3  # how fractionally close the SECOND-TO-LAST beta grid point comes to betamax_scalar
    map = {}
    #array: idx to tuple of e_val, array
    for e_launch in e_vals:
        map[round(e_launch, 5)] = {}
        for cosf_launch in cosf:
            betamax = (1. - e_launch ** 2) / (2. * (1. + e_launch * cosf_launch))
            betamax = stabfac * betamax
            dNdbeta = lambda beta: beta ** betapow * (1 - beta) ** 1.5 * (
                    1 - e_launch ** 2 - 2 * beta * (1 + e_launch * cosf_launch)) ** -1.5

            epsmax = np.log10(betamax - betamin)
            epsmin = np.log10(how_close * betamax)
            eps = np.linspace(epsmax, epsmin, n_beta_grid - 1)
            beta = betamax - 10. ** eps
            beta = np.append(beta, betamax)

            norm = sint.quad(dNdbeta, betamin, betamax)[0]
            Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0] / norm
            Nbeta = np.vectorize(Nbeta)
            invNbeta = si.interp1d(Nbeta(beta), beta)
            map[round(e_launch,5)][round(cosf_launch,5)] = invNbeta
    # print()
    # print(e_vals)
    # print(cosf)
    return e_vals, f_vals, map

def OrbTimeCorr_MidOptimized_Background(num_e=10, max_e=0.2, num_f=10, n_beta_grid=50, betapow=1.5, betamin=0.001, stabfac=0.997):

    #return inverse CDF for e, cosf combo
    #hashmap: e_val to list of cosf values

    e_vals = np.linspace(0, max_e, num_e)
    cosf = np.linspace(-1, 1, num_f)
    f_vals = np.arccos(cosf)
    for i in range(len(f_vals)):
        if nr.uniform(0, 1) > 0.5:
            f_vals[i] = 2 * np.pi - f_vals[i]
    cosf = np.cos(f_vals)
    how_close = 1.e-3  # how fractionally close the SECOND-TO-LAST beta grid point comes to betamax_scalar
    map = {}
    #array: idx to tuple of e_val, array
    for e_launch in e_vals:
        map[round(e_launch, 5)] = {}
        for cosf_launch in cosf:
            betamax = (1. - e_launch ** 2) / (2. * (1. + e_launch * cosf_launch))
            betamax = stabfac * betamax
            dNdbeta = lambda beta: beta ** betapow * (1 - beta) ** 1.5 * (
                    1 - e_launch ** 2 - 2 * beta * (1 + e_launch * cosf_launch)) ** -1.5

            epsmax = np.log10(betamax - betamin)
            epsmin = np.log10(how_close * betamax)
            eps = np.linspace(epsmax, epsmin, n_beta_grid - 1)
            beta = betamax - 10. ** eps
            beta = np.append(beta, betamax)

            norm = sint.quad(dNdbeta, betamin, betamax)[0]
            Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0] / norm
            Nbeta = np.vectorize(Nbeta)
            invNbeta = si.interp1d(Nbeta(beta), beta)
            map[round(e_launch,5)][round(cosf_launch,5)] = invNbeta
    # print()
    # print(e_vals)
    # print(cosf)
    return e_vals, f_vals, map

def OrbTimeCorr_Background_ApplyBeta(invNbeta, a_launch, e_launch, I_launch, Omega_launch, omega_launch, cosf_launch, sinf_launch,
                             beta_per_launch, n_beta_grid=50, betapow=1.5, betamin=0.001, stabfac=0.997):
    
    for j in range(len(e_launch)): # loop over different directions
        e_scalar = e_launch[j]  # store scalars for slight speed advantage
        cosf_scalar = cosf_launch[j]
        betamax_scalar = betamax[j]

        # dNdbeta should account for non-zero e and non-unity cosf_scalar
        dNdbeta = lambda beta: beta ** betapow * (1 - beta) ** 1.5 * (
                    1 - e_scalar ** 2 - 2 * beta * (1 + e_scalar * cosf_scalar)) ** -1.5

        # create logarithmic grid, focused near betamax
        how_close = 1.e-3  # how fractionally close the SECOND-TO-LAST beta grid point comes to betamax_scalar
        epsmax = np.log10(betamax_scalar - betamin)
        epsmin = np.log10(how_close * betamax_scalar)
        eps = np.linspace(epsmax, epsmin, n_beta_grid - 1)
        beta = betamax_scalar - 10. ** eps
        beta = np.append(beta, betamax_scalar)

        norm = sint.quad(dNdbeta, betamin, betamax_scalar)[0]
        Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0] / norm
        Nbeta = np.vectorize(Nbeta)

        invNbeta = si.interp1d(Nbeta(beta), beta)
        beta_segment = invNbeta(nr.uniform(0, 1, beta_per_launch))
        beg = j * beta_per_launch
        end = beg + beta_per_launch

        a_dust[beg:end] = (1 - beta_segment) * a_launch[j] * (1 - e_launch[j] ** 2) / (
                    1 - e_launch[j] ** 2 - 2 * beta_segment * (1 + e_launch[j] * cosf_launch[j]))
        omega_dust[beg:end] = omega_launch[j] + np.arctan2(beta_segment * sinf_launch[j],
                                                           e_launch[j] + beta_segment * cosf_launch[j])
        e_dust[beg:end] = np.sqrt(
            e_launch[j] ** 2 + 2 * beta_segment * e_launch[j] * cosf_launch[j] + beta_segment ** 2) / (1 - beta_segment)
        I_dust[beg:end] = I_launch[j]
        Omega_dust[beg:end] = Omega_launch[j]
        beta_dust[beg:end] = beta_segment
    return a_dust, e_dust, I_dust, Omega_dust, omega_dust, beta_dust 

def OrbTimeCorr_MidOptimized(a_launch, e_launch, I_launch, Omega_launch, omega_launch, cosf_launch, sinf_launch,
                             beta_per_launch, n_beta_grid=50, betapow=1.5, betamin=0.001, stabfac=0.997, beta_limit=1):
    betamax = (1. - e_launch ** 2) / (2. * (1. + e_launch * cosf_launch))
    for i in range(len(betamax)):
        if betamax[i] > beta_limit:
            betamax[i] = beta_limit
    betamax = stabfac * betamax
    Nlaunch = len(a_launch)
    total_dust_particles = Nlaunch * beta_per_launch
    a_dust = np.zeros(total_dust_particles)
    e_dust = np.zeros(total_dust_particles)
    I_dust = np.zeros(total_dust_particles)
    Omega_dust = np.zeros(total_dust_particles)
    omega_dust = np.zeros(total_dust_particles)
    beta_dust = np.zeros(total_dust_particles)
    for j in range(len(e_launch)): # loop over different directions
        e_scalar = e_launch[j]  # store scalars for slight speed advantage
        cosf_scalar = cosf_launch[j]
        betamax_scalar = betamax[j]

        # dNdbeta should account for non-zero e and non-unity cosf_scalar
        dNdbeta = lambda beta: beta ** betapow * (1 - beta) ** 1.5 * (
                    1 - e_scalar ** 2 - 2 * beta * (1 + e_scalar * cosf_scalar)) ** -1.5

        # create logarithmic grid, focused near betamax
        how_close = 1.e-3  # how fractionally close the SECOND-TO-LAST beta grid point comes to betamax_scalar
        epsmax = np.log10(betamax_scalar - betamin)
        epsmin = np.log10(how_close * betamax_scalar)
        eps = np.linspace(epsmax, epsmin, n_beta_grid - 1)
        beta = betamax_scalar - 10. ** eps
        beta = np.append(beta, betamax_scalar)

        norm = sint.quad(dNdbeta, betamin, betamax_scalar)[0]
        Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0] / norm
        Nbeta = np.vectorize(Nbeta)

        invNbeta = si.interp1d(Nbeta(beta), beta)
        beta_segment = invNbeta(nr.uniform(0, 1, beta_per_launch))
        beg = j * beta_per_launch
        end = beg + beta_per_launch

        a_dust[beg:end] = (1 - beta_segment) * a_launch[j] * (1 - e_launch[j] ** 2) / (
                    1 - e_launch[j] ** 2 - 2 * beta_segment * (1 + e_launch[j] * cosf_launch[j]))
        omega_dust[beg:end] = omega_launch[j] + np.arctan2(beta_segment * sinf_launch[j],
                                                           e_launch[j] + beta_segment * cosf_launch[j])
        e_dust[beg:end] = np.sqrt(
            e_launch[j] ** 2 + 2 * beta_segment * e_launch[j] * cosf_launch[j] + beta_segment ** 2) / (1 - beta_segment)
        I_dust[beg:end] = I_launch[j]
        Omega_dust[beg:end] = Omega_launch[j]
        beta_dust[beg:end] = beta_segment
    return a_dust, e_dust, I_dust, Omega_dust, omega_dust, beta_dust




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


def get_inverse_CDF_stab(e=0, cosf=1, betapow=1.5, betamin=0.001, betamax=1, stabfac=0.997, precision=100):
    dNdbeta = lambda beta: beta ** betapow * (1 - beta) ** 1.5 * (1 - e ** 2 - 2 * beta * (1 + e * cosf)) ** -1.5
    f_betamax = (1 - e ** 2) / 2. / (1 + e * cosf) * stabfac
    if f_betamax > betamax:
        f_betamax = stabfac * betamax

    norm = sint.quad(dNdbeta, betamin, f_betamax)[0]
    Nbeta = lambda beta: sint.quad(dNdbeta, betamin, beta)[0] / norm
    Nbeta = np.vectorize(Nbeta)
    invNbeta = si.interp1d(Nbeta(np.linspace(betamin * 0.8, f_betamax, precision)),
                           np.linspace(betamin * 0.8, f_betamax, precision))
    return invNbeta


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

