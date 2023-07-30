import numpy as np
import numpy.random as nr
import scipy.integrate as sint
import scipy.interpolate as si


# Samples betas. Only samples a limited number due to speed constraints
def OrbTimeCorr(a_launch, e_launch, I_launch, Omega_launch, omega_launch, cosf_launch, sinf_launch,
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


def OrbTimeCorrBackground(e_set, num_f=10, n_beta_grid=50, betapow=1.5, betamin=0.001, stabfac=0.997):
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
    return e_vals, f_vals, map





