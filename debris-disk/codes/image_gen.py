from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from numpy import power,sin,cos
from plt_macro import mac_setup as setup
from matplotlib.patches import Ellipse

X, Y = 400,400
RES = 0.1
Xgrid, Ygrid = np.meshgrid(np.arange(-X/2,X/2,RES), np.arange(-Y/2,Y/2,RES))
def pol2cart(r,theta):
    return np.array( [r*cos(theta),r*sin(theta) ] ).T

def radius(_thetas,_as,_es,_oms):
    return _as * (1 - power(_es, 2))/(1+_es*cos(_thetas-_oms))

def getmaski(ad, ed, alt, az):
    y_ell = ad*np.sin(alt*np.pi/180.)
    ycent = ed*ad*np.sin(alt*np.pi/180.)*np.sin(az*np.pi/180.)
    x_ell = ad*np.sqrt(1-ed**2)
    xcent = ed*ad*np.sin((az+90)*np.pi/180.)
    
    Yscl = (Ygrid+ycent)/y_ell
    Xscl = (Xgrid+xcent)/x_ell

    return np.where( np.power(Yscl,2.)+np.power(Xscl,2.) < 1.)

setup()
thetas = np.linspace(0,2*np.pi,1000)

r1s = radius(thetas,33,0.7,0)
r2s = radius(thetas,33*1.33,0.556,0)

x1s,y1s = pol2cart(r1s,thetas).T
x2s,y2s = pol2cart(r2s,thetas).T

X,Y = 400,400

alts = [-3] #np.ceil(np.linspace(-2.5,-87.5,16)) #range(-2.5,-95,-5)
i=1
for alt in alts:
    az = 270
    print(alt)
    fstr = "single"
    #fstr = 'moth_0I_uniformM_evo'
    fstr2 = 'moth_0I'

    r1s = radius(thetas,33,0.7,0)
    r2s = radius(thetas,33*1.33,0.556,0)

    x1s,y1s = pol2cart(r1s,thetas).T
    x2s,y2s = pol2cart(r2s,thetas).T

    fig,ax = plt.subplots(1,1,sharey=True,figsize = (12,12))
    Im = np.loadtxt("C:/Users/jon88/Desktop/debris-disk/images/imgrid/%ix%i/%s_Alt%04i_Az%04i.txt"%(X/2,Y/2,fstr,alt,az))
    a_d, e_d, I_d, Om_d, om_d = np.loadtxt("C:/Users/jon88/Desktop/debris-disk/dustorbit/%s_dustorbit.txt"%fstr2, unpack=True, usecols=(0,1,2,3,4))


    #ai = min(a_d)
    #ei = e_d[np.argmin(a_d)]
    #bi = ai*np.sqrt(1-ei**2 )*np.cos(np.pi/2-alt*np.pi/180)

    #e = Ellipse(xy = (0,0), width=2*ai, height = 2*bi,zorder=100)

    #ax.add_artist(e)
    #e.set_facecolor('pink')

    Ims = ndimage.filters.gaussian_filter(Im, 2, mode='nearest')


    maski = getmaski(a_d[np.argmin(a_d)], e_d[np.argmin(a_d)], alt, az)
    Ims[maski] = np.nan

    cm = plt.cm.hot
    cm.set_bad('black',1.)

    ax.imshow(np.sqrt(Ims), cmap=cm, origin='lower', \
            vmin=0,vmax=4, extent=[-X / 2, X / 2, -Y / 2, Y / 2])


    #ax.contour(np.sqrt(Ims),levels=[0.1,0.2,0.3,0.4,0.5,0.6] ,colors=['1'],extent=[-X / 2, X / 2, -Y / 2, Y / 2])

    #ax.plot(x1s,y1s,color='0',ls='dashed',lw=2.5)
    #ax.plot(x2s,y2s,color='0',ls='dashed',lw=2.5)

    ax.set_xlabel('x (AU)')
    ax.set_ylabel('y (AU)')

    ax.set_aspect('equal')
    plt.tight_layout()
    #fig.subplots_adjust(wspace=0)
    #fig.savefig('alt'+str(i))
    i+=2
    print(i)
    plt.show()


'''
setup()
thetas = np.linspace(0,2*np.pi,1000)

r1s = radius(thetas,33,0.7,0)
r2s = radius(thetas,33*1.33,0.556,0)

x1s,y1s = pol2cart(r1s,thetas).T
x2s,y2s = pol2cart(r2s,thetas).T

X,Y = 400,400
alt = 10
az = 90
fstr0 = 'moth_0I'
fstr1 = 'moth_0I_e07_evo'

fig,(ax0,ax1) = plt.subplots(1,2,sharey=True,figsize = (8,4))
Im0 = np.loadtxt("C:/Users/jon88/Desktop/debris-disk/images/imgrid/%ix%i/%s_Alt%04i_Az%04i.txt"%(X/2,Y/2,fstr0,alt,az))
Im1 = np.loadtxt("C:/Users/jon88/Desktop/debris-disk/images/imgrid/%ix%i/%s_Alt%04i_Az%04i.txt"%(X/2,Y/2,fstr1,alt,az))
Ims0 = ndimage.filters.gaussian_filter(Im0, 2, mode='nearest')
Ims1 = ndimage.filters.gaussian_filter(Im1, 2, mode='nearest')

ax0.imshow(np.sqrt(Ims0), cmap=plt.cm.hot, origin='lower', \
           vmin=0,vmax=80, extent=[-X / 2, X / 2, -Y / 2, Y / 2])

ax1.imshow(np.sqrt(Ims1), cmap=plt.cm.hot, origin='lower', \
           vmin=0,vmax=80, extent=[-X / 2, X / 2, -Y / 2, Y / 2])

ax0.set_xlabel('x (AU)')
ax1.set_xlabel('x (AU)')
ax0.set_ylabel('y (AU)')

ax0.set_aspect('equal')
ax1.set_aspect('equal')
plt.tight_layout()
#fig.subplots_adjust(wspace=0)
plt.show()
'''
'''
lvls = np.logspace(-3,1.5,100)

fix, ax = plt.subplots()
Image = np.loadtxt("C:/Users/jon88/Desktop/debris-disk/images/imgrid/%ix%i/%s_Alt%04i_Az%04i.txt"%(X/2,Y/2,fstr,alt,az))
Imsmooth = ndimage.filters.gaussian_filter(Image, 2, mode='nearest')

print(np.max(np.sqrt(Imsmooth)))

ax.imshow(np.sqrt(Imsmooth), cmap=plt.cm.hot, origin='lower', \
           vmin=0,vmax=110, extent=[-X / 2, X / 2, -Y / 2, Y / 2])
#ax.contourf(np.sqrt(Imsmooth), cmap=plt.cm.hot, origin='lower', levels=lvls, vmin=0,\
#           extent=[-X / 2, X / 2, -Y / 2, Y / 2])

#ax.plot(x1s,y1s,color='1',ls='dashed')
#ax.plot(x2s,y2s,color='1',ls='dashed')

plt.axis('equal')
plt.show()
'''
#apsidal precession frequency
#gauss's equations to calculate orbital changes after kicks