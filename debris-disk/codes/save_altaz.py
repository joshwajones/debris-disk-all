"""
Save scattered light image
Array of alt+azimuth

@author: Eve J. Lee
Mar. 31st 2016
"""

import numpy as np
import matplotlib.pyplot as pl
import scatter_image
import sys


maxa = np.float(sys.argv[1])
d = np.float(sys.argv[2]) 
fstr = sys.argv[3]
alt = float(sys.argv[4])
az = float(sys.argv[5])
ar = np.float(sys.argv[6])
if len(sys.argv) > 7:
    #print(len(sys.argv))
    #print(sys.argv)
    Nd = int(sys.argv[7])
    #print(Nd)
else:
    Nd = 100
verbose=False
every_print=0
if len(sys.argv) > 8:
    every_print = int(sys.argv[8])
    verbose = True

include_depth = False
include_azfirst = False
alt_only = True
if include_depth:
    Image_alt, Image_az, col, opt = scatter_image.MakeImage("../dustorbit/%s_dustorbit.txt" % fstr, aspect_ratio=ar,
                                                            resolution=0.1, obsincl=alt, maxa=maxa, d=d, obsazim=az,
                                                            Ndust=Nd,
                                                            depth=include_depth, fixbeta=0, verbose=verbose,
                                                            every_x_print=every_print, include_azfirst=include_azfirst)

    np.savetxt("../images/imgrid/%ix%i/%s_Alt%04i_Az%04i_azfirst.txt" % (maxa, maxa * ar, fstr, alt, az), Image_az)
    np.savetxt("../images/imgrid/%ix%i/%s_Alt%04i_Az%04i_coldensity.txt" % (maxa, maxa * ar, fstr, alt, az), col)
    np.savetxt("../images/imgrid/%ix%i/%s_Alt%04i_Az%04i_optdepth.txt" % (maxa, maxa * ar, fstr, alt, az), opt)
elif alt_only:
    Image_alt = scatter_image.MakeImage_altonly("../dustorbit/%s_dustorbit.txt"%fstr, aspect_ratio=ar,
                                                            resolution=0.1, obsincl=alt,maxa=maxa,d=d,obsazim=az,Ndust=Nd,
                                                            fixbeta=0, verbose=verbose, every_x_print=every_print,)

else:
    Image_alt, Image_az = scatter_image.MakeImage("../dustorbit/%s_dustorbit.txt" % fstr, aspect_ratio=ar,
                                                            resolution=0.1, obsincl=alt, maxa=maxa, d=d, obsazim=az,
                                                            Ndust=Nd,
                                                            include_depth=include_depth, fixbeta=0, verbose=verbose,
                                                            every_x_print=every_print)

    np.savetxt("../images/imgrid/%ix%i/%s_Alt%04i_Az%04i_azfirst.txt" % (maxa, maxa * ar, fstr, alt, az), Image_az)
np.savetxt("../images/imgrid/%ix%i/%s_Alt%04i_Az%04i.txt"%(maxa, maxa*ar,fstr,alt,az), Image_alt)

