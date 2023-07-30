import numpy as np
import scatter_image
import sys
import time
import re

# typical command line: python save_altaz.py 400 10 name alt az 1 100 10000 4 1000
start_time = time.time()

if len(sys.argv) == 2: # reading from file
    filename = "../image_gen/" + sys.argv[1] + ".txt"
    read_data = {}
    with open(filename, 'r') as fileobject:
        for line in fileobject:
            label, data = line.split('=')
            read_data[label.strip()] = data.strip()
    maxa = np.float(read_data["maxa"])
    d = np.float(read_data["d"])
    fstr = read_data["fstr"]
    if "alt_az_pairs" in read_data:
        altaz_pairs = read_data["alt_az_pairs"].strip()
        altaz_pairs = re.split("\)[\s ,;]*\(", altaz_pairs)
        for i, pair in enumerate(altaz_pairs):
            pair = pair.strip("()")
            pair = re.split(",\s*|\s*,", pair)
            altaz_pairs[i] = [float(num) for num in pair]
    else:
        alt = float(read_data["alt"])
        az = float(read_data["az"])
        altaz_pairs = [[alt, az]]

    ar = np.float(read_data["ar"])
    Nd = int(read_data.get("Nd", 100))
    every_print = int(read_data.get("print_every_x_dust"), 0)
    verbose = (every_print > 0)
    SPF = int(read_data.get("SPF", 4))
    wavelength = float(read_data.get("wavelength", 100))
else: # from command line
    maxa = np.float(sys.argv[1])
    d = np.float(sys.argv[2])
    fstr = sys.argv[3]
    alt = float(sys.argv[4])
    az = float(sys.argv[5])
    ar = np.float(sys.argv[6])
    if len(sys.argv) > 7:
        Nd = int(sys.argv[7])
    else:
        Nd = 100
    verbose=False
    every_print=0
    if len(sys.argv) > 8:
        every_print = int(sys.argv[8])
        verbose = True
    SPF = 4
    if len(sys.argv) > 9:
        SPF = int(sys.argv[9])
    wavelength = 1000
    if len(sys.argv) > 10:
        wavelength = int(sys.argv[10])

include_depth = False
include_azfirst = False
alt_only = True
for alt, az in altaz_pairs:
    if include_depth: # mostly deprecated. Saves alt-first image, az-first image, column density.
        Image_alt, Image_az, col, opt = scatter_image.MakeImage("../dustorbit/%s_dustorbit.txt" % fstr, aspect_ratio=ar,
                                                                resolution=0.1, obsincl=alt, maxa=maxa, d=d, obsazim=az,
                                                                Ndust=Nd,
                                                                depth=include_depth, fixbeta=0, verbose=verbose,
                                                                every_x_print=every_print, include_azfirst=include_azfirst)

        np.savetxt("../images/imgrid/%ix%i/%s_Alt%04i_Az%04i_azfirst.txt" % (maxa, maxa * ar, fstr, alt, az), Image_az)
        np.savetxt("../images/imgrid/%ix%i/%s_Alt%04i_Az%04i_coldensity.txt" % (maxa, maxa * ar, fstr, alt, az), col)
        np.savetxt("../images/imgrid/%ix%i/%s_Alt%04i_Az%04i_optdepth.txt" % (maxa, maxa * ar, fstr, alt, az), opt)
    elif alt_only: # What I used almost exclusively and what I recommend. Only the alt-first image.
        Image_alt = scatter_image.MakeImage_altonly("../dustorbit/%s_dustorbit.txt"%fstr, aspect_ratio=ar,
                                                                resolution=0.1, obsincl=alt,maxa=maxa,d=d,obsazim=az,Ndust=Nd,
                                                                fixbeta=0, verbose=verbose, every_x_print=every_print, SPF=SPF, wavelength=wavelength)

    else:
        Image_alt, Image_az = scatter_image.MakeImage("../dustorbit/%s_dustorbit.txt" % fstr, aspect_ratio=ar,
                                                                resolution=0.1, obsincl=alt, maxa=maxa, d=d, obsazim=az,
                                                                Ndust=Nd,
                                                                include_depth=include_depth, fixbeta=0, verbose=verbose,
                                                                every_x_print=every_print)

        np.savetxt("../images/imgrid/%ix%i/%s_Alt%04i_Az%04i_azfirst.txt" % (maxa, maxa * ar, fstr, alt, az), Image_az)
    np.savetxt("../images/imgrid/%ix%i/%s_Alt%04i_Az%04i.txt"%(maxa, maxa*ar,fstr,alt,az), Image_alt)
    print("\n")
end_time = time.time()
total_time = end_time - start_time
print("Total time:   ", total_time)
