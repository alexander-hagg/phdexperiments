"""
frenet.py

Spectral analysis of frenet variables given an x,y input. 

"""
import sys
sys.path.insert(1, '/home/ryi/projects_py/lib')

import numpy as np
from scipy import signal

# Custom.
import interparc as ia
from twodimnav import curvature
from psd import psdfull

# Functions for data transformation. 

def calc_dnds(x,y):
    """ dnds calculates the derivative dN/ds, with respect to the normal vector. 
    """
    # Deal with arrays again.
    x = np.asarray(x)
    y = np.asarray(y)
    x = x.reshape(len(x),)
    y = y.reshape(len(y),)
    # Calculate curvature.
    kappa = curvature(x,y)
    tanvec = np.diff(x) + np.diff(y)*1j 
    return -kappa*(tanvec[0:-1] + tanvec[1:])/2 

def calc_dtds(x,y):
    """ dtds calculates the derivative dT/ds, with respect to the tangent.
    """
    # Deal with arrays again.
    x = np.asarray(x)
    y = np.asarray(y)
    x = x.reshape(len(x),)
    y = y.reshape(len(y),)
    # Calculate curvature.
    kappa = curvature(x,y)
    normvec = -np.diff(y) + np.diff(x)*1j
    return kappa*(normvec[0:-1] + normvec[1:])/2 

# Functions for spectral analysis.

def frenet_periodogram_i(x, y, l=0):
    """ periodogram calculates per. for inputs x, y, after interpolating over l=len(x) (default).
    """
    if l==0:
        l = len(x)
    interp_coords = ia.interparc(l, x, y)
    x_i = interp_coords[:,0]
    y_i = interp_coords[:,1]
    # Calculate frenet. 
    dtds_i = calc_dtds(x_i, y_i)
    steps_i = np.sqrt(np.diff(x_i, axis=0)**2 + np.diff(y_i, axis=0)**2)[:-1]
    fs = 1./steps_i[0]
    # Calculate periodogram.
#    per_f, per_p = psdfull(dtds_i, fs)
    per_f, per_p = signal.periodogram(dtds_i, fs, return_onesided=False)
    per_f = np.fft.fftshift(per_f)
    per_p = np.fft.fftshift(per_p)
    return per_f, per_p, fs, x_i, y_i
