""" 
psd.py

This contains a single function, psdfull, which calculates the full power spectral density of an input series.  
"""

import numpy as np

def psdfull(x, Fs, window='hann'):
    """ 
    psdfull

    Calculates the discrete power spectral density per unit time (we almost always want per unit time!).
    'x' is the input series. 'Fs' is sampling frequency. 
    """
    # Rename variables to conventional names. 
    # Length of the input series.
    N = len(x)
    tmax = N/Fs
   
    # Calculate the FFT. By default use the hanning window. 
    if window.lower == 'hann':
        xdft = np.fft.fft(x*np.hanning(len(x))) 
    else:
        xdft = np.fft.fft(x)
    # Shift the fft to be symmetric.
    xdft = np.fft.fftshift(xdft)
    # Calculate the power spectrum. 
    psdx = abs(xdft)**2

    # Manually create an fftshifted index vector (this could be done more efficiently, but whatever). 
    freq_index= np.fft.fftshift(np.arange(N))
    # Find the zero frequency index.
    zero_index = np.where(freq_index == 0)[0][0]
    # Let all frequencies left of this be negative frequencies.
    leftside = -np.arange(zero_index)-1
    # Insert this vector (flipped) into the freq_index vector.
    freq_index[0:zero_index] = leftside[::-1]
   
    # Convert to frequency vector.
    freq = freq_index/tmax

    # Rescale. See page 389 of 'Numerical Methods'. First the default dft sum ignores the sampling frequency and just sums the terms.
    psdx = psdx/(Fs**2)
    # We generally want the signal to be independent of the number of samples (psd per unit time). So divide by tmax.
    psdx = psdx/tmax

    return freq, psdx
