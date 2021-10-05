#!/usr/bin/env python
# coding: utf-8

NUMBER_OF_MODES = 6


import numpy as np
import matplotlib.pyplot as plt
from pydmd import DMD, CDMD

snapshots=[]

for i in range(30000,100000,500):
    u0 = np.load('fig{}.npy'.format(i))[0]
    u1 = np.load('fig{}.npy'.format(i))[1]
    snapshots+=[np.sqrt(u0*u0+u1*u1)]

dmd = CDMD(NUMBER_OF_MODES)
dmd.fit(snapshots)

modes = dmd.modes
for i in range(NUMBER_OF_MODES):
    plt.imshow(modes[:,i].reshape(300,200).real.T,cmap='seismic')
    plt.savefig('mode{}.pdf'.format(i))


