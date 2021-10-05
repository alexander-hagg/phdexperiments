"""
anglevslength.py

A module that contains functions to calculate angle as a function of arclength.

"""
import numpy as np
import twodimnav as tdn

def anglevslength(x, y, n):
    """ anglevslength(x, y)
    Calculates the angle as a function of length for coordinates x and y over n points.
    """
    # Calculate arclength of x and y.
    arc = tdn.arclength(x, y)
    # Average 'arc' so it is the same size as the angles.
    for i in range(n-1):
        arc = (arc[0:-1] + arc[1:])/2.
    # Initialize vector for angle storage.
    angles = np.zeros(len(x))
    # Calculate the angle between successive segments. 
    # Deal with edge effects by iteratively looping and re-inserting. This is inefficient, but n, x, and y are usually small so who cares?  
    dx = x[n-1:] - x[0:-n+1]
    dy = y[n-1:] - y[0:-n+1]
    angles = np.arctan(dy/dx)
    return arc, angles
