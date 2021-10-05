"""
frenet.py

Allows for easy calculation of frenet frame variables given an x,y input. 

"""
import sys
sys.path.insert(1, '/home/ryi/projects_py/lib')
import numpy as np
from twodimnav import curvature

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
