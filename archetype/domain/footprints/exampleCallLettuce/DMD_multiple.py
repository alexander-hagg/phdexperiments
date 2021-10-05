#!/usr/bin/env python
# coding: utf-8

NUMBER_OF_MODES = 24


import numpy as np
import matplotlib.pyplot as plt
from pydmd import DMD, CDMD
import scipy.io

######################################################
# High Frequency Long Run
snapshots=[]

for i in range(50000,100000,50):
    u0 = np.load('fig{}.npy'.format(i))[0]
    u1 = np.load('fig{}.npy'.format(i))[1]
    snapshots+=[np.sqrt(u0*u0+u1*u1)]

dmd = CDMD(NUMBER_OF_MODES)
dmd.fit(snapshots)

modes = dmd.modes
scipy.io.savemat('modes_long_high.mat',{'m':modes[:,0].reshape(300,200).real.T})

plt.imshow(modes[:,0].reshape(300,200).real.T,cmap='seismic')
plt.savefig('mode_long_high_freq.pdf')


######################################################
# High Frequency 
snapshots=[]

for i in range(90000,100000,50):
    u0 = np.load('fig{}.npy'.format(i))[0]
    u1 = np.load('fig{}.npy'.format(i))[1]
    snapshots+=[np.sqrt(u0*u0+u1*u1)]

dmd = CDMD(NUMBER_OF_MODES)
dmd.fit(snapshots)

modes = dmd.modes
scipy.io.savemat('modes_high.mat',{'m':modes[:,0].reshape(300,200).real.T})

plt.imshow(modes[:,0].reshape(300,200).real.T,cmap='seismic')
plt.savefig('mode_high_freq.pdf')

######################################################
# Medium Frequency 
snapshots=[]

NUMBER_OF_MODES = 12

for i in range(70000,100000,500):
    u0 = np.load('fig{}.npy'.format(i))[0]
    u1 = np.load('fig{}.npy'.format(i))[1]
    snapshots+=[np.sqrt(u0*u0+u1*u1)]

dmd = CDMD(NUMBER_OF_MODES)
dmd.fit(snapshots)
modes = dmd.modes
scipy.io.savemat('modes_medium.mat',{'m':modes[:,0].reshape(300,200).real.T})

plt.imshow(modes[:,0].reshape(300,200).real.T,cmap='seismic')
plt.savefig('mode_medium_freq.pdf')

######################################################
# Low Frequency 

snapshots=[]

NUMBER_OF_MODES = 6

for i in range(50000,100000,2000):
    u0 = np.load('fig{}.npy'.format(i))[0]
    u1 = np.load('fig{}.npy'.format(i))[1]
    snapshots+=[np.sqrt(u0*u0+u1*u1)]

dmd = CDMD(NUMBER_OF_MODES)
dmd.fit(snapshots)
modes = dmd.modes
scipy.io.savemat('modes_low.mat',{'m':modes[:,0].reshape(300,200).real.T})

plt.imshow(modes[:,0].reshape(300,200).real.T,cmap='seismic')
plt.savefig('mode_low_freq.pdf')


