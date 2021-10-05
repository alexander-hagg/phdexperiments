"""
twodimnav.py

Functions that help determine distances, closest points, whether
points are on lines between other points in two dimensions. 

"""

import numpy as np
import stdize 

# Functions for linear geometry. 
def arclength(x, y):
    """ Calculates the arclength along coordinates x and y. x and y must be numpy arrays. 
    """
    x = stdize.stdvec(x)
    y = stdize.stdvec(y)
    # Calculate the distances for each point. 
    arc = np.hypot(np.diff(x), np.diff(y))
    arc = np.append(0, arc)
    arc = np.cumsum(arc, dtype=float)
    return arc
def curvature(x, y, n=0):
    """ Calculate the curvature for x, y coords. 
    """
    # Make sure the format is correct.
    x = stdize.stdvec(x)
    y = stdize.stdvec(y)
    # Calculate first derivative.
    d1 = np.sqrt(np.diff(y)**2 + np.diff(x)**2)
    dx = np.diff(x)/d1
    dy = np.diff(y)/d1
    d1x = (dx[0:-1] + dx[1:])/2
    d1y = (dy[0:-1] + dy[1:])/2
    # Calculate second derivative.
    d2 = (d1[0:-1] + d1[1:])/2
    d2x = np.diff(dx)/d2
    d2y = np.diff(dy)/d2  
    # Smooth the derivatives (optional). 
    if n:
        d1x = smooth(d1x, n)
        d1y = smooth(d1y, n)
        d2x = smooth(d2x, n)
        d2y = smooth(d2y, n)
    # Calculate curvature.
    return (d1x*d2y - d1y*d2x)/((d1x**2 + d1y**2)**(3/2)) 
def distance(a,b):
    """ Calculate distance between two points, two lists of points, or a list and a point. This requires 2d numpy arrays as input, even for two points. 
    """
    # Resize if inputs are vectors.
    if len(a.shape) == 1:
        a = a.reshape(1, a.shape[0])
    if len(b.shape) == 1:
        b = b.reshape(1, b.shape[0])    
    # Calculate distance.
    dist = np.sqrt(np.subtract(a[:,0], b[:,0])**2 + np.subtract(a[:,1], b[:,1])**2)
    return dist
def find_closest_points(x, points, n):
    """ Find the 'n' closest points to 'x' in 'points'. 
    """
    d = distance(points, x)
    ds = np.argsort(d)
    sorted_points = points.take(ds, axis=0)
    return sorted_points[0:n]
def is_between(a,c,b):
    """ Define function for determining if point is on a vertice. 
    """
    tol = 1E-8
    return -tol < distance(a,c) + distance(c,b) - distance(a,b) < tol
def on_branch(x, branch):
    """ Module to determine if boundary condition is on a given branch. 
    """
    # Find closest vertices to 'x'. 
    vertices = find_closest_points(x, branch, 2)
    # Check if 'x' is in between these vertices.
    return is_between(vertices[0], x, vertices[1])
def stepsize(x, y):
    """ Calculates the stepsizes between two vectors x and y.
    """
    x = stdize.stdvec(x)
    y = stdize.stdvec(y)
    steps = np.hypot(np.diff(x), np.diff(y))
    return steps

# Functions for mesh editing. 
def symmetric_point(x, y, distance):
    """ Returns two points given a vertex: one on either side of vertex a distance 'distance' away in a direction perpendicular to vertex. 
    """
    # Calculate the angle of the segment.
    dy = y[-1]-y[0]
    dx = x[-1]-x[0]
    angle = np.arctan2(dy, dx)
    # Calculate coordinates of point +- 90 degrees offset. 
    pl = np.array([x[1] + distance*np.cos(angle - np.pi/2.), y[1] + distance*np.sin(angle - np.pi/2.)])
    pr = np.array([x[1] + distance*np.cos(angle + np.pi/2.), y[1] + distance*np.sin(angle + np.pi/2.)])
    return pl, pr
def symmetric_list(branch, distance):
    # Module to calculate the symmetric points for an entire branch of points.
    points_left = np.array([]).reshape(0,2)
    points_right = np.array([]).reshape(0,2)
    # Loop through 'branch' and add points on either side.
    for i, point in enumerate(branch):
        if i>0:
            x = branch[i-1:min(i+1, len(branch)),0]
            y = branch[i-1:min(i+1, len(branch)),1]
            pl, pr = symmetric_point(x, y, distance)
            points_left = np.append(points_left, [pl], axis=0)
            points_right = np.append(points_right, [pr], axis=0)
    return points_left, points_right

# Fitting functions.
def rma(x, y):
    """ rma(x,y): reduced major axis regression on x and y."""
    Sxy = np.cov(x,y)
    Sxy = Sxy[0,1]

    rho = Sxy/(np.std(x, ddof=1)*np.std(y, ddof=1))
    b = np.std(y, ddof=1)/np.std(x, ddof=1)
    beta = np.sign(rho)*b
    alpha = np.mean(y) - beta*np.mean(x)
    return beta, alpha, rho

# Smoothing function.
def smooth(x, n):
    """ smooth(x,n): smooths x over a moving window of size n. """
    x = stdize.stdvec(x) 
    return np.convolve(x, np.ones((n,))/n, mode='full')
    
