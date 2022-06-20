"""
Check a single scattered light image

@author: Eve J. Lee
Apr. 21st 2016
"""

import numpy as np
import matplotlib.pyplot as pl
from scipy import ndimage
import read_disk_input
import plt_macro

#X, Y = 800, 800
X, Y  = 400, 400
Xgrid, Ygrid = np.meshgrid(np.arange(-X/2,X/2), np.arange(-Y/2,Y/2))

def getmaski(ad, ed, alt, az):
    y_ell = ad*np.sin(alt*np.pi/180.)
    ycent = ed*ad*np.sin(alt*np.pi/180.)*np.sin(az*np.pi/180.)
    x_ell = ad*np.sqrt(1-ed**2)
    xcent = ed*ad*np.sin((az+90)*np.pi/180.)

    Yscl = (Ygrid+ycent)/y_ell
    Xscl = (Xgrid+xcent)/x_ell

    return np.where(Yscl**2+Xscl**2 < 1)

#alt, az = 90, 0
alt, az = 45, 45
pl.clf()
#fstr = "single_inner_planet_e0.7_Ifree0_efree0_betadistrb1.5_bmin0.001_bmax1.0_orbcorr_morelaunch"
fstr = "moth_big"
#fstr = "test_2"
Image = np.loadtxt("../images/imgrid/%ix%i/%s_Alt%04i_Az%04i.txt"%(X/2,Y/2,fstr,alt,az))
Imsmooth = ndimage.filters.gaussian_filter(Image, 2, mode='nearest')

a_p, e_p, I_p, Om_p, om_p = np.loadtxt("../parentorbit/%s_orbit.txt"%fstr, unpack=True)
a_d, e_d, I_d, Om_d, om_d = np.loadtxt("../dustorbit/%s_dustorbit.txt"%fstr, unpack=True, usecols=(0,1,2,3,4))

pl.imshow(np.sqrt(Imsmooth), cmap=pl.cm.hot, origin='lower', \
          vmin=0, vmax=20, extent=[-X/2, X/2, -Y/2, Y/2])

maski = getmaski(a_d[np.argmin(a_d)], e_d[np.argmin(a_d)], alt, az)
Imsmooth[maski] = np.nan

#pl.contour(Xgrid, Ygrid, Imsmooth, levels=[1.6], colors='white', origin='lower')
pl.plot(0,0,'yo',markersize=5)
#pl.xlim((-150,150))
pl.ylim((-300,300))
pl.show()

