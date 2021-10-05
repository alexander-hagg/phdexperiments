"""
stdize.py

Operations on objects that convert them to a standard workable format. 

"""
import numpy as np

def stdvec(x):
    """ stdvec returns a 1D numpy array, given a list or 2d numpy array.
    """
    x = np.asarray(x)
    x = x.reshape(len(x),)
    return x
