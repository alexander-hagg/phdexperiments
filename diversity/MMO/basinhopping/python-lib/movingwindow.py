import numpy as np
from numpy import convolve

def movingwindow(values, window):
""" Averages 'values' over a moving window of size 'window' (uses numpy 'convolve').
Idea came from 'gordoncluster.wordpress.com'.

Returns a vector of the size of value - window + 1. I.e. edges are ignored. Edges should be added back in.  
"""
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma
