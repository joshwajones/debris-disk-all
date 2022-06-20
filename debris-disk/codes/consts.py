"""
Program to load constants
all in cgs

@author: Eve Lee
Mar. 6th 2013
"""

import numpy as np

#fundamental constants
G = 6.67e-8 #gravitational constant
h = 6.63e-27 #planck's constant
h_ev = 4.14e-15 #planck's constant in eV.s
hbar = h/(2.*np.pi)
k = 1.38e-16 #Boltzmann's constant
sb = 5.67e-5 #Stefan-Boltzmann constant
e = 4.80e-10 #electron charge in statcoulombs
e_coul = 1.6e-19 #electron charge in coulombs
mH = 1.67e-24 #mass of proton
m_e = 9.11e-28 #mass of electron
c = 3e10 #speed of light
sig_T = 6.65e-25 #Thomson scattering cross section in cgs
a0 = 5.29e-9 #Bohr radius in cm
arad = 4.*sb/c #radiation constant

#Diatomic hydrogen molecule
Hbondlen = 74e-10
H2momI = 0.5*mH*Hbondlen**2.
H2freq = 1.32e14
H2bondE = 4.5167 #in eV
H2rotE = hbar**2./2./H2momI
H2vibE = hbar*H2freq

#conversion
s2yr = 3600.*24.*365.25 #actually yr to sec
s2Myr = s2yr*1e6 #actually Myr to sec
pc2cm = 3e18
kpc2cm = 3e21
au2cm = 1.5e13
km2cm = 1e5
jy2cgs = 1e-23
as2rad = np.pi/180./3600.
eV2ergs = 1.6e-12
barns2cgs = 1e-24
rad2arcsec = 206264.8

#solar properties
Msun = 2e33
Lsun = 4e33
Rsun = 7e10
Teff_sun = 5800.
Tc_sun = 1.44e7
rhoc_sun = 89.9
rho_avg_sun = 1.
sol_X = 0.75
sol_Y = 0.23
sol_Z = 0.02
sol_mu = 0.6

#Earth properties
Mearth = 5.97e27
Rearth = 6371.*1e5

#Jupiter properties
Mjup = 1.9e30
Rjup = 6.99e9

#cosmic abundance
cosmo_X = 0.7
cosmo_Y = 0.28
cosmo_Z = 0.02
cosmo_mu = 0.62

#cosmology
H0 = 70. #Hubble's constant in km/s/Mpc

#Radiation
Q2M = 6.34e46 #ionizing Q to mass in Msun
alpha_H = 3.57e-13 #cm^-3 s^-1, H recombination coefficient
Tion = 7000. #K temperature of ionized gas
Otmf = 4. #Myr lifetime of O stars in clusters

#GMC Mstar probes
g_Mstar2Lbol = 1.879093e+03
g_M2Q = 6.263861e+46
g_tmsLbol = 11.82843 #in Myr
g_tmsQ = 3.866884 #in Myr
