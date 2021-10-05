"""
curvature.py 

Functions that calculate curvature of a line. 
"""

import numpy as np

def curvature(x, y, n):
    # Calculates the curvature over n points for a moving window in x and y.
    # Calculate the first derivative.
    d = np.sqrt(np.diff(x)**2. + np.diff(y)**2.)
    dx = np.diff(x)/d
    dy = np.diff(y)/d

    # Calculate the second derivative.
    d1 = (d[:-1] + d[1:])/2.
    d2x = np.diff(dx)/d1
    d2y = np.diff(dy)/d1

    # Adjust the first derivative to have the same number of points.
    d1x = (dx[:-1] + dx[1:])/2.
    d1y = (dy[:-1] + dy[1:])/2.

    # Calculate curvature. 
    kappa = (d1x*d2y - d1y*d2x)/((d1x**2. + d1y**2.)**(3./2.))
    return kappa
