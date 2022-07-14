"""
Function to make a single 
scattered light image

@author: Eve J. Lee
Feb. 24th 2016
"""
import numpy as np
import prob_position_orbit as ppo
import consts
import time

HG = lambda mu, g: (1./4./np.pi)*(1-g**2)/(1+g**2-2*g*mu)**1.5
compHG = lambda mu: 0.643*HG(mu,0.995) + 0.176*HG(mu,0.665) + 0.181*HG(mu,0.035)    #G Ring fit

use_compHG = True


def MakeImage_altonly(dustfile, d=10, maxa=100., aspect_ratio=1., resolution=0.05, L=1., g=0.5, Ndust=100, obsincl=5.,
              obsazim=0., fixbeta=0., verbose=False, every_x_print=1):
    # dustfile = name of the file containint dust orbital parameters
    # d = distance to the disk in pc
    # maxa = maximum horizontal lengthscale covered by the image in AU
    # aspect_ratio = height/width
    # resolution = pixelscale of the image in arcsec
    # L = luminosity of the host star in solar luminosity
    # g = HG scattering phase function anisotropy parameter [-1, 1] (-1 backward scattering, 0 isotropic, 1 forward scattering)
    # Ndust = number of dust particles per orbit
    # obsincl, obsazim = viewing angle in degrees

    # viewing angle
    obsincl *= np.pi / 180.
    obsazim *= np.pi / 180.

    # set up image grid
    Npix_x = int((2 * maxa / d) / resolution)
    Npix_y = int((2 * maxa / d) * aspect_ratio / resolution)
    imx = np.linspace(-maxa / d, maxa / d, Npix_x)
    imy = np.linspace(-maxa / d * aspect_ratio, maxa / d * aspect_ratio, Npix_y)
    image_alt = np.zeros((Npix_y, Npix_x))

    # load dust orbits
    if fixbeta == 0:
        a, e, I, Omega, omega, beta = np.loadtxt(dustfile, unpack=True)
    else:
        a, e, I, Omega, omega = np.loadtxt(dustfile, unpack=True)
        beta = np.ones(len(a)) * fixbeta

    for i in range(len(a)):  # for each dust grain
        if verbose and i % every_x_print == 0:
            print("Inc:%4.2f, Az:%4.2f, %i out of %i" % (obsincl * 180. / np.pi, obsazim * 180. / np.pi, i + 1, len(a)))

        # samples Npts = Ndust = 100 (default) true anomalies by selecting mean anomalies uniformly
        f, fweight = ppo.OutputPosition(e[i], Npts=Ndust)
        # radius at each point
        Rd = (a[i] * (1 - e[i] ** 2) / (1 + e[i] * np.cos(f)))

        # xyz coords of drawn points in 3d space
        # on unit sphere: not scaled by R
        Xdraw = np.cos(Omega[i]) * np.cos(omega[i] + f) - np.sin(Omega[i]) * np.sin(omega[i] + f) * np.cos(I[i])
        Ydraw = np.sin(Omega[i]) * np.cos(omega[i] + f) + np.cos(Omega[i]) * np.sin(omega[i] + f) * np.cos(I[i])
        Zdraw = np.sin(omega[i] + f) * np.sin(I[i])

        Xd_alt = Xdraw * np.cos(obsazim) - Ydraw * np.sin(obsazim)
        Yd_alt = Xdraw * np.sin(obsazim) * np.cos(obsincl) + Ydraw * np.cos(obsazim) * np.cos(obsincl) - Zdraw * np.sin(
            obsincl)
        Zd_alt = Xdraw * np.sin(obsazim) * np.sin(obsincl) + Ydraw * np.cos(obsazim) * np.sin(obsincl) + Zdraw * np.cos(
            obsincl)

        if not use_compHG:
            intensity_alt = (L * consts.Lsun / 4. / np.pi / (Rd * consts.au2cm) ** 2) * HG(-Yd_alt, g) * (
                        beta[i] * 10) ** -2
        else:
            intensity_alt = (L * consts.Lsun / 4. / np.pi / (Rd * consts.au2cm) ** 2) * compHG(-Yd_alt) * (
                        beta[i] * 10) ** -2
        # for every sample point associated with an orbit
        for r, x_alt, z_alt, idust_alt in zip(Rd, Xd_alt, Zd_alt, intensity_alt):
            Xind_alt = np.searchsorted(imx, x_alt * r / d)
            Zind_alt = np.searchsorted(imy, z_alt * r / d)
            if Zind_alt < Npix_y and Xind_alt < Npix_x:
                image_alt[Zind_alt, Xind_alt] += idust_alt
    return image_alt

def MakeImage(dustfile, d=10, maxa=100., aspect_ratio=1., resolution=0.05, L=1., g=0.5, Ndust=100, obsincl=5., obsazim=0., fixbeta=0., include_depth=False, verbose=False, every_x_print=0):
    #dustfile = name of the file containint dust orbital parameters
    #d = distance to the disk in pc
    #maxa = maximum horizontal lengthscale covered by the image in AU
    #aspect_ratio = height/width
    #resolution = pixelscale of the image in arcsec
    #L = luminosity of the host star in solar luminosity
    #g = HG scattering phase function anisotropy parameter [-1, 1] (-1 backward scattering, 0 isotropic, 1 forward scattering)
    #Ndust = number of dust particles per orbit
    #obsincl, obsazim = viewing angle in degrees

    #viewing angle
    obsincl *= np.pi/180.
    obsazim *= np.pi/180.
    
    #set up image grid
    Npix_x = int((2*maxa/d)/resolution)
    Npix_y = int((2*maxa/d)*aspect_ratio/resolution)
    imx = np.linspace(-maxa/d, maxa/d, Npix_x)
    imy = np.linspace(-maxa/d*aspect_ratio, maxa/d*aspect_ratio, Npix_y)
    image_alt = np.zeros((Npix_y, Npix_x))
    if include_depth:
        col_density = np.zeros((Npix_y, Npix_x))
        opt_depth = np.zeros((Npix_y, Npix_x))
    image_az = np.zeros((Npix_y, Npix_x))

    #load dust orbits
    if fixbeta == 0:
        a, e, I, Omega, omega, beta = np.loadtxt(dustfile, unpack=True)    
    else:
        a, e, I, Omega, omega = np.loadtxt(dustfile, unpack=True)    
        beta = np.ones(len(a))*fixbeta
 
    for i in range(len(a)): #for each dust grain
        if verbose and i % every_x_print == 0: 
                print("Inc:%4.2f, Az:%4.2f, %i out of %i"%(obsincl*180./np.pi, obsazim*180./np.pi, i+1, len(a)))

        #samples Npts = Ndust = 100 (default) true anomalies by selecting mean anomalies uniformly
        f, fweight = ppo.OutputPosition(e[i], Npts=Ndust)
        #radius at each point
        Rd = (a[i]*(1-e[i]**2)/(1+e[i]*np.cos(f)))

        #xyz coords of drawn points in 3d space
        #on unit sphere: not scaled by R
        Xdraw = np.cos(Omega[i])*np.cos(omega[i]+f)-np.sin(Omega[i])*np.sin(omega[i]+f)*np.cos(I[i])
        Ydraw = np.sin(Omega[i])*np.cos(omega[i]+f)+np.cos(Omega[i])*np.sin(omega[i]+f)*np.cos(I[i])
        Zdraw = np.sin(omega[i]+f)*np.sin(I[i])

        Xd_alt = Xdraw*np.cos(obsazim)-Ydraw*np.sin(obsazim)
        Yd_alt = Xdraw*np.sin(obsazim)*np.cos(obsincl)+Ydraw*np.cos(obsazim)*np.cos(obsincl)-Zdraw*np.sin(obsincl)
        Zd_alt = Xdraw*np.sin(obsazim)*np.sin(obsincl)+Ydraw*np.cos(obsazim)*np.sin(obsincl)+Zdraw*np.cos(obsincl)
        
        Xd_az = Xdraw*np.cos(obsazim)-Ydraw*np.cos(obsincl)*np.sin(obsazim)+Zdraw*np.sin(obsincl)*np.sin(obsazim)
        Yd_az = Xdraw*np.sin(obsazim)+Ydraw*np.cos(obsincl)*np.cos(obsazim)-Zdraw*np.sin(obsincl)*np.cos(obsazim)
        Zd_az = Ydraw*np.sin(obsincl)+Zdraw*np.cos(obsincl)

        if not use_compHG:
            intensity_alt = (L*consts.Lsun/4./np.pi/(Rd*consts.au2cm)**2)*HG(-Yd_alt, g)*(beta[i]*10)**-2
            
            intensity_az = (L*consts.Lsun/4./np.pi/(Rd*consts.au2cm)**2)*HG(-Yd_az, g)*(beta[i]*10)**-2
        else:
            intensity_alt = (L*consts.Lsun/4./np.pi/(Rd*consts.au2cm)**2)*compHG(-Yd_alt)*(beta[i]*10)**-2
            intensity_az = (L*consts.Lsun/4./np.pi/(Rd*consts.au2cm)**2)*compHG(-Yd_az)*(beta[i]*10)**-2

        #for every sample point associated with an orbit
        for r, x_alt, z_alt, idust_alt, x_az, z_az, idust_az in zip(Rd, Xd_alt, Zd_alt, intensity_alt, Xd_az, Zd_az, intensity_az):
            #print(r, x_alt, z_alt, idust_alt)
            # if i==0:
            #     print("r:     ", r)
            #     print("x:     ", x_alt)
            #     print("z:     ", z_alt)
            Xind_alt = np.searchsorted(imx, x_alt*r/d)
            Zind_alt = np.searchsorted(imy, z_alt*r/d)
            if Zind_alt < Npix_y and Xind_alt < Npix_x:
                image_alt[Zind_alt, Xind_alt] += idust_alt
                if include_depth:
                    col_density[Zind_alt, Xind_alt] += 1
                    opt_depth[Zind_alt, Xind_alt] += 1 / (beta[i] ** 2)
            
            Xind_az = np.searchsorted(imx, x_az*r/d)
            Zind_az = np.searchsorted(imy, z_az*r/d)
            if Zind_az < Npix_y and Xind_az < Npix_x:
                image_az[Zind_az, Xind_az] += idust_az
    if include_depth: 
        return image_alt, image_az, col_density, opt_depth
    else:
        return image_alt, image_az
    
