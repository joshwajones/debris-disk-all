"""
Calculate forced mutual inclination
for secular perturbation by two planets

@author: Eve J. Lee
Dec. 21st 2017
"""

import numpy as np
import matplotlib.pyplot as pl
import matplotlib.patches as patches
import fn_freq as ff
import consts

adisk = np.logspace(0, 2, 1000)
B = np.zeros(len(adisk))
B1 = np.zeros(len(adisk))
B2 = np.zeros(len(adisk))

Mps = np.array([0.78, 1])*consts.Mjup/consts.Mearth
aps = np.array([3.48, 48])

for i in range(len(adisk)):
    B[i] = ff.B(adisk[i], Mps, aps, Mstar=0.82*consts.Msun)
    B1[i] = ff.Bj(adisk[i], Mps[0], aps[0], Mstar=0.82 * consts.Msun)
    B2[i] = ff.Bj(adisk[i], Mps[1], aps[1], Mstar=0.82 * consts.Msun)

Bmat = ff.Bmatrix(aps[0], aps[1], Mps[0], Mps[1], Mstar=0.82*consts.Msun)
eigval, eigvec = np.linalg.eig(Bmat)
eigvec = np.transpose(eigvec)

EpsErib_i = 89
EpsErib_imin = 89-42
EpsErib_imax = 89+42

EpsEri_iout = 33.

init_I1 = 0.
init_I2 = 20. #45 #EpsErib_i-EpsEri_iout

mu_1 = B1*init_I1*eigvec[0][0]+B2*init_I2*eigvec[1][0]
mu_2 = B1*init_I1*eigvec[0][1]+B2*init_I2*eigvec[1][1]

age = 800 #Myrs
p0 = -(mu_1/(B-eigval[0]))*np.sin(eigval[0]*age*consts.s2Myr) - \
     (mu_2/(B-eigval[1]))*np.sin(eigval[1]*age*consts.s2Myr)
q0 = -(mu_1/(B-eigval[0]))*np.cos(eigval[0]*age*consts.s2Myr) - \
     (mu_2/(B-eigval[1]))*np.cos(eigval[1]*age*consts.s2Myr)

I = np.sqrt(p0**2+q0**2)

fig, ax = pl.subplots(1, 1)
ax.plot(adisk, I, 'k-', lw=2)
ax.axhline(EpsErib_i-EpsEri_iout, lw=2)
ax.axhspan(EpsErib_imin-EpsEri_iout, EpsErib_imax-EpsEri_iout, color='orange', alpha=0.5)
#ax.add_patch(patches.Rectangle((62.6, EpsErib_imin-EpsEri_iout), \
#                               75.9-62.6, EpsErib_imax-EpsErib_imin, \
#                               facecolor='red', alpha=0.5))
ax.axvspan(62.6, 75.9, color='blue', alpha=0.2)

pl.tight_layout()
pl.show()


