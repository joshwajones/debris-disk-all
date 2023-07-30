import math
import pdb

import numpy as np
import matplotlib.pyplot as pl
from scipy import ndimage
import time
import astropy
from astropy.io import fits
from astropy.visualization import astropy_mpl_style


# For each filename in filenames, saves each alt-az az pair as a FITS file (if FITS) and displays (if show)
def save_files(filenames, FITS=True, txt=True, show=False, alts=[], azes=[], smooth=False, maxa=400):
    X = Y = maxa
    Xgrid, Ygrid = np.meshgrid(np.arange(-X/2,X/2), np.arange(-Y/2,Y/2))

    imgrid_path = "../images/imgrid/" + str(X) + "x" + str(Y) + "/"
    for filename in filenames:
        for alt in alts:
            for az in azes:
                try:
                    alt_angle = get_angle(alt)
                    az_angle = get_angle(az)
                    name = imgrid_path + filename + "_Alt" + alt_angle + f"_Az{az_angle}.txt"
                    Image = np.loadtxt(name)
                    Imsmooth = ndimage.filters.gaussian_filter(Image, 2, mode='nearest')
                    if smooth:
                        Image = Imsmooth
                    if txt:
                        display = np.power(Imsmooth, 1. / 2.)
                        pl.figure(1)
                        pl.clf()
                        pl.imshow(display, cmap=pl.cm.hot, origin='lower', \
                                  vmin=np.amin(display), vmax=20, extent=[-X / 2, X / 2, -Y / 2, Y / 2])
                        pl.colorbar(label="Pixel value", orientation="horizontal", fraction=0.05, pad=0.1)
                        pl.plot(0, 0, 'yo', markersize=5)
                        pl.ylim((-300, 300))
                        if show:
                            pl.pause(1)
                    if FITS:
                        hdu = fits.PrimaryHDU(Image)
                        fullFile = filename + "_Alt" + alt_angle + f"_Az{az_angle}"
                        if smooth:
                            fullFile = fullFile + "_smooth"
                        save_file_name = "../FITS/" + fullFile
                        hdu.writeto(save_file_name + ".fits", overwrite=True)
                except:
                    print("FAILED ON " + str(alt) + ", " + str(az))


# Converts int angle to string format
def get_angle(angle):
    angle = int(angle)
    suffix = str((abs(angle)))
    sign = "0"
    if angle < 0:
        sign = "-"
    second_char = ""
    if abs(angle) < 100:
        second_char = "0"

    third_char = ""
    if abs(angle) < 10:
        third_char = "0"

    return sign + second_char + third_char + suffix


def animation(filename="", wait_time=3, alts=[], azes=[]):
    X, Y = 800, 800
    Xgrid, Ygrid = np.meshgrid(np.arange(-X/2,X/2), np.arange(-Y/2,Y/2))

    imgrid_path = "..images/imgrid/"
    for alt in alts:
        for az in azes:
            try:
                alt_angle = get_angle(alt)
                az_angle = get_angle(az)
                name = imgrid_path + "400x400/" + filename + "_Alt" + alt_angle + "_Az" + az_angle + ".txt"
                Image = np.loadtxt(name)
                Imsmooth = ndimage.filters.gaussian_filter(Image, 2, mode='nearest')

                display = np.power(Imsmooth, 1. / 2.)
                pl.figure(1)
                pl.clf()
                pl.imshow(display, cmap=pl.cm.hot, origin='lower', \
                          vmin=np.amin(display), vmax=20, extent=[-X / 2, X / 2, -Y / 2, Y / 2])
                pl.colorbar(label="Pixel value", orientation="horizontal", fraction=0.05, pad=0.1)
                pl.title(f"{alt_angle}, {az_angle}l]")
                pl.plot(0, 0, 'yo', markersize=5)
                pl.ylim((-300, 300))
                pl.pause(wait_time)
            except:
                pass

def display_img(filename, alt=90, az=90, size=400, ar=1., vmax=20, sqrted=False, smooth=True, imgrid_path = "../images/imgrid/",
                xL=-400, xR=400, ylim=400):
    """
    Check a single scattered light image

    @author: Eve J. Lee
    Apr. 21st 2016
    """
    X, Y = 2*size, 2*ar * size
    Xgrid, Ygrid = np.meshgrid(np.arange(-X/2,X/2), np.arange(-Y/2,Y/2))
    alt = get_angle(alt)
    az = get_angle(az)
    def getmaski(ad, ed, alt, az):
        y_ell = ad*np.sin(alt*np.pi/180.)
        ycent = ed*ad*np.sin(alt*np.pi/180.)*np.sin(az*np.pi/180.)
        x_ell = ad*np.sqrt(1-ed**2)
        xcent = ed*ad*np.sin((az+90)*np.pi/180.)

        Yscl = (Ygrid+ycent)/y_ell
        Xscl = (Xgrid+xcent)/x_ell

        return np.where(Yscl**2+Xscl**2 < 1)

    pl.clf()
    ysize = int(size * ar)
    size = str(size)
    ysize = str(ysize)
    imgrid_path = imgrid_path + size +"x" + ysize + "/"
    fullFile = filename + "_Alt" + alt + "_Az" + az
    name = imgrid_path + fullFile + ".txt"
    Image = np.loadtxt(name)
    Imsmooth = ndimage.filters.gaussian_filter(Image, 2, mode='nearest')
    if smooth:
        display = Imsmooth
    else:
        display = Image
    if sqrted:
        display = np.power(display, 1./2.)
    pl.imshow(display, cmap=pl.cm.hot, origin='lower', \
              vmin=0, vmax=vmax, extent=[-X/2, X/2, -Y/2, Y/2], aspect='auto')
    #pl.colorbar(label="Pixel value", orientation="horizontal", fraction=0.05, pad=0.1)
    pl.title(f"{alt}, {az}")
    pl.plot(0,0,'yo',markersize=5)
    pl.ylim((-ylim,ylim))
    pl.xlim((xL, xR))
    pl.show()

# Used to overlay multiple images
def combine_imgs(file1, file2, alts, azes_1, azes_2, saveName=""):
    X = Y= 800
    imgrid_path = "../images/imgrid/400x400/"
    if not saveName:
        saveName = file1 + file2
    for alt in alts:
        alt_angle = get_angle(alt)
        for i in range(len(azes_1)):
            az1_angle = get_angle(azes_1[i])
            az2_angle = get_angle(azes_2[i])
            try:
                name1 = imgrid_path + file1 + "_Alt" + alt_angle + f"_Az{az1_angle}.txt"
                name2 = imgrid_path + file2 + "_Alt" + alt_angle + f"_Az{az2_angle}.txt"

                Im1 = np.loadtxt(name1)
                Im2 = np.loadtxt(name2)
                Image = Im1 + Im2
                fullName = f"../images/imgrid/400x400/{saveName}_Alt{alt_angle}_Az{az1_angle}.txt"
                np.savetxt(fullName, Image)
            except:
                print("FAILED ON: alt = " + str(alt) + " with az1 = " + str(az1_angle) + " and az2 = " + str(az2_angle))
    save_files(filenames=[saveName], FITS=True, txt=True, show=False, alts=alts, azes=azes_1)


# PUT ACTUAL ACTION YOU WANT TO TAKE HERE
# EXAMPLE OF SAVING FILES:
filenames = ["thermal_high_ecc_test"]
alts = [10]
azes = [11]
# save_files(filenames=filenames, alts=alts, azes=azes)


# EXAMPLE OF COMBINING TWO IMAGES:
# combine_imgs("back_coll_3_thermal", "inc_coll_8_thermal", [-3], [150], [150], "together_thermal_20-20-1")

# EXAMPLE OF DISPLAYING:
display_img("back_coll_3_thermal", alt=-3, az=150)

