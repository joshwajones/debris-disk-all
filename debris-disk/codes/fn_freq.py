"""
Functions to calculate 
precession frequencies

@author: Eve J. Lee
July 21st 2015
"""
import numpy as np
import consts
import scipy.integrate as sint

def alpha(a, ap):
    return np.where(ap < a, ap/a, a/ap)

def malpha(a, ap):
    return np.where(ap < a, 1, a/ap)

def b(s, j, alpha):
    integrand = lambda psi, j, s, alpha: np.cos(j * psi) / (1 - 2 * alpha * np.cos(psi) + alpha ** 2.) ** s
    try:
        bfinal = np.zeros(len(alpha))
        for i in range(len(alpha)):
            bfinal[i] = (1. / np.pi) * sint.quad(integrand, 0, 2 * np.pi, args=(j, s, alpha[i]))[0]
    except TypeError:
        bfinal = (1. / np.pi) * sint.quad(integrand, 0, 2 * np.pi, args=(j, s, alpha))[0]

    return bfinal

def A(a, Mp, ap, Mstar=consts.Msun):
    n = np.sqrt(consts.G*Mstar/(a*consts.au2cm)**3.)
    if type(Mp) == type(np.array([])):
        val1 = 0.25*n*(Mp[0]*consts.Mearth/Mstar)*alpha(a,ap[0])*malpha(a,ap[0])*b(1.5, 1, alpha(a,ap[0]))
        val2 = 0.25*n*(Mp[1]*consts.Mearth/Mstar)*alpha(a,ap[1])*malpha(a,ap[1])*b(1.5, 1, alpha(a,ap[1]))
        return val1+val2
    else:
        return 0.25*n*(Mp*consts.Mearth/Mstar)*alpha(a,ap)*malpha(a,ap)*b(1.5, 1, alpha(a,ap))
      
def B(a, Mp, ap, Mstar=consts.Msun):
    n = np.sqrt(consts.G*Mstar/(a*consts.au2cm)**3.)
    if type(Mp) == type(np.array([])):
        val1 = -0.25*n*(Mp[0]*consts.Mearth/Mstar)*alpha(a,ap[0])*malpha(a,ap[0])*b(1.5, 1, alpha(a,ap[0]))
        val2 = -0.25*n*(Mp[1]*consts.Mearth/Mstar)*alpha(a,ap[1])*malpha(a,ap[1])*b(1.5, 1, alpha(a,ap[1]))
        return val1+val2
    else:
        return -0.25*n*(Mp*consts.Mearth/Mstar)*alpha(a,ap)*malpha(a,ap)*b(1.5, 1, alpha(a,ap))

def Aj(a, Mp, ap, Mstar=consts.Msun):
    n = np.sqrt(consts.G*Mstar/(a*consts.au2cm)**3.)
    return -0.25*n*(Mp*consts.Mearth/Mstar)*alpha(a,ap)*malpha(a,ap)*b(1.5,2,alpha(a,ap))

def Bj(a, Mp, ap, Mstar=consts.Msun):
    n = np.sqrt(consts.G*Mstar/(a*consts.au2cm)**3.)
    return 0.25*n*(Mp*consts.Mearth/Mstar)*alpha(a,ap)*malpha(a,ap)*b(1.5,1,alpha(a,ap))

def Amatrix(a1, a2, M1, M2, Mstar=consts.Msun):
    n1 = np.sqrt(consts.G*Mstar/(a1*consts.au2cm)**3.)
    n2 = np.sqrt(consts.G*Mstar/(a2*consts.au2cm)**3.)
    
    A11 = 0.25*n1*(M2*consts.Mearth/(Mstar+M1*consts.Mearth))*alpha(a1,a2)*malpha(a1,a2)*b(1.5,1,alpha(a1,a2))
    A12 = -0.25*n1*(M2*consts.Mearth/(Mstar+M1*consts.Mearth))*alpha(a1,a2)*malpha(a1,a2)*b(1.5,2,alpha(a1,a2))
    A21 = -0.25*n2*(M1*consts.Mearth/(Mstar+M2*consts.Mearth))*alpha(a2,a1)*malpha(a2,a1)*b(1.5,2,alpha(a2,a1))
    A22 = 0.25*n2*(M1*consts.Mearth/(Mstar+M2*consts.Mearth))*alpha(a2,a1)*malpha(a2,a1)*b(1.5,1,alpha(a2,a1))

    return np.array([[A11,A12],[A21,A22]])

def Bmatrix(a1, a2, M1, M2, Mstar=consts.Msun):
    n1 = np.sqrt(consts.G*Mstar/(a1*consts.au2cm)**3.)
    n2 = np.sqrt(consts.G*Mstar/(a2*consts.au2cm)**3.)
    
    B11 = -0.25*n1*(M2*consts.Mearth/(Mstar+M1*consts.Mearth))*alpha(a1,a2)*malpha(a1,a2)*b(1.5,1,alpha(a1,a2))
    B12 = 0.25*n1*(M2*consts.Mearth/(Mstar+M1*consts.Mearth))*alpha(a1,a2)*malpha(a1,a2)*b(1.5,1,alpha(a1,a2))
    B21 = 0.25*n2*(M1*consts.Mearth/(Mstar+M2*consts.Mearth))*alpha(a2,a1)*malpha(a2,a1)*b(1.5,1,alpha(a2,a1))
    B22 = -0.25*n2*(M1*consts.Mearth/(Mstar+M2*consts.Mearth))*alpha(a2,a1)*malpha(a2,a1)*b(1.5,1,alpha(a2,a1))

    return np.array([[B11,B12],[B21,B22]])





