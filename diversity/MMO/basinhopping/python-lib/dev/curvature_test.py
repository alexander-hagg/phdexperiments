""" 
curvature_test.py

Testing script for 'curvature.py' package. 
"""

import curvature as crv
import numpy as np

x = np.array([0,1,2,3,4,5])
y = x**2
n = 1

a = crv.curvature(x, y, n)

print(a)
