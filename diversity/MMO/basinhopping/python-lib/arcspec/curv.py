""" 
curv.py


"""
import sys
sys.path.insert(1, '/home/ryi/projects_py/lib')

from astropy.stats import LombScargle
import numpy as np
from scipy import signal

# Custom. 
import interparc as ia
from psd import psdfull
from twodimnav import curvature

def periodogram_i(x, y, l=0):
    """ periodogram_i calculates per. for inputs x, y, after interpolating over l=len(x) (default).
    """
    if l==0:
        l = len(x)
    interp_coords = ia.interparc(l, x, y)
    x_i = interp_coords[:,0]
    y_i = interp_coords[:,1]
    # Calculate curvature. 
    curv = curvature(x_i, y_i)
    steps = np.sqrt(np.diff(x_i, axis=0)**2 + np.diff(y_i, axis=0)**2)[:-1]
    arc = np.cumsum(steps)
    fs = 1./steps[0]
    # Calculate periodogram.
    per_f, per_p = psdfull(curv, fs)
    return per_f, per_p, fs 

def lombs(x, y):
    """ lombs calculates LombScargle for inputs x, y.
    """
    # Calculate curvature. 
    curv = curvature(x, y)
    steps = np.sqrt(np.diff(x, axis=0)**2 + np.diff(y, axis=0)**2)[:-1]
    arc = np.cumsum(steps)
    # Calculate LS.
    ls_f, ls_p = LombScargle(arc, curv).autopower()
    return ls_f, ls_p

def arccurv(x, y):
    """ Just calculates the curvature and the arclengths corresponding to those curvatures. 
    """
    curv = curvature(x, y)
    steps = np.sqrt(np.diff(x, axis=0)**2 + np.diff(y, axis=0)**2)[:-1]
    arc = np.cumsum(steps)
    return arc, curv

def arccurv_i(x, y, l=0):
    """ Calculates the curvature and associated arclength for an interpolated channel.
    """
    if l==0:
        l = len(x)
    interp_coords = ia.interparc(l, x, y)
    x_i = interp_coords[:,0]
    y_i = interp_coords[:,1]
    # Calculate curvature. 
    curv = curvature(x_i, y_i)
    steps = np.sqrt(np.diff(x_i, axis=0)**2 + np.diff(y_i, axis=0)**2)[:-1]
    arc = np.cumsum(steps)
    return arc, curv
