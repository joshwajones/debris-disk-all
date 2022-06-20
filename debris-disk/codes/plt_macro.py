"""
Plot macro

@author: Eve Lee
Mar. 6th 2013
"""

from matplotlib import rc
from matplotlib import rcParams

def mac_setup(fontsize=20, usetex=True, xtick=20, xtmaj=15, xtmin=10, xtmaj_w=2, xtmin_w=2, ytick=20, ytmaj=15, ytmin=10, ytmaj_w=2, ytmin_w=2, axes=20, legend=18, weight='bold', boldTeX=True):
## Set up plot parameters
    rc('font', size=fontsize, weight=weight, **{'family':'serif'})
    rc('text', usetex=usetex)
    if boldTeX:
        rcParams['text.latex.preamble']=[r'\boldmath']
    rc('xtick', labelsize=xtick)
    rc('xtick.major', size=xtmaj, width=xtmaj_w)
    rc('xtick.minor', size=xtmin, width=xtmin_w)
    rc('ytick', labelsize=ytick)
    rc('ytick.major', size=ytmaj, width=ytmaj_w)
    rc('ytick.minor', size=ytmin, width=ytmin_w)
    rc('axes', labelsize=axes)
    rc('legend', fontsize=legend)


