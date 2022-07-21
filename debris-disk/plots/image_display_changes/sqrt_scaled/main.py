import math
import numpy as np
import matplotlib.pyplot as pl
from scipy import ndimage


#all important image generation happens in smoothed_img()


# get total SB at radius from center
def get_tot(im):
    vals = np.zeros(400)
    y_center = len(im) // 2
    x_center = len(im[0]) // 2
    for y in range(len(im)):
        for x in range(len(im)):
            delta_y = y - y_center
            delta_x = x - x_center
            r = np.sqrt(delta_x**2 + delta_y**2)
            if r<400:
                vals[math.floor(r)] += im[y][x]
    return vals

# used to fit a plot (optical depth, column density) to a line
def plot_line(filename, lbl, smoothed=True, pl_line=False):
    im = np.loadtxt(filename)
    if smoothed:
        im = ndimage.filters.gaussian_filter(im, 2, mode='nearest')
    row = im[len(im) // 2]
    yvals = np.log10(row[(len(row) + 1) // 2:])

    xvals = [i for i in range(1, len(yvals) + 1)]
    xvals = np.log10(xvals)
    pl.plot(xvals, yvals, label=lbl)
    if pl_line:
        x_line_vals = xvals[80:]
        b = 5.75 + .1
        m = -1.5
        y_line_vals = [b + m * x for x in x_line_vals]
        pl.plot(x_line_vals, y_line_vals, label=f"slope={m}")


# show the raw, unfiltered image
def raw_img():
    Image = np.loadtxt("/Users/sjosh/PycharmProjects/Research/img_2/debris-disk/images/imgrid/400x400/background_14_Alt0000_Az0000.txt")
    pl.imshow(Image)
    pl.show()


def smoothed_img():
    """
    Check a single scattered light image

    @author: Eve J. Lee
    Apr. 21st 2016
    """
    X, Y = 800, 800
    Xgrid, Ygrid = np.meshgrid(np.arange(-X/2,X/2), np.arange(-Y/2,Y/2))

    def getmaski(ad, ed, alt, az):
        y_ell = ad*np.sin(alt*np.pi/180.)
        ycent = ed*ad*np.sin(alt*np.pi/180.)*np.sin(az*np.pi/180.)
        x_ell = ad*np.sqrt(1-ed**2)
        xcent = ed*ad*np.sin((az+90)*np.pi/180.)

        Yscl = (Ygrid+ycent)/y_ell
        Xscl = (Xgrid+xcent)/x_ell

        return np.where(Yscl**2+Xscl**2 < 1)

    pl.clf()
    imgrid_path = "/Users/sjosh/PycharmProjects/Research/img_2/debris-disk/images/imgrid/"
    #Image = np.loadtxt(imgrid_path + "400x400/moth_repro_10mil_Alt-010_Az0090.txt")
    #Image = np.loadtxt(imgrid_path + "400x400/moth_repro_10mil_Alt0090_Az0090.txt")
    Image = np.loadtxt(imgrid_path + "400x400/moth_repro_10mil_Alt-003_Az0090.txt")
    Imsmooth = ndimage.filters.gaussian_filter(Image, 2, mode='nearest')

    display = np.power(Imsmooth, 1./2.)
    pl.imshow(display, cmap=pl.cm.hot, origin='lower', \
              vmin=np.amin(display), vmax=np.amax(display), extent=[-X/2, X/2, -Y/2, Y/2])
    pl.colorbar(label="Pixel value", orientation="horizontal", fraction=0.05, pad=0.1)

    #pl.contour(Xgrid, Ygrid, Imsmooth, levels=[1.6], colors='white', origin='lower')
    pl.plot(0,0,'yo',markersize=5)
    #pl.xlim((-150,150))
    pl.ylim((-300,300))
    pl.show()

smoothed_img()
