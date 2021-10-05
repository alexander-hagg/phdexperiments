"""
bin_data.py 

Different routines for binning data.

scipy.stats.binned_statistic does the same thing as bin_x (but they require edges as input, not centers). Ugh what a waste of time. 

"""

from cleanarrays import standardarray
import numpy as np

def bin_x(x, y, gx):
    """ bin_x bins the input data 'x' into centers at 'gx' and calculates the mean of 'y' in those bins. 
    """
    # Make sure the format is correct.
    x = standardarray(x)
    y = standardarray(y)
    gx = standardarray(gx)

    # Remove NaNs here.
    good_x = np.where(~np.isnan(x))[0]
    good_y = np.where(~np.isnan(y))[0]
    good_idx = np.intersect1d(good_x, good_y)
    x = x[good_idx]
    y = y[good_idx]

    # Deal with edge factors.
    binwidth = np.diff(gx)
    gx = np.hstack((gx[0] - binwidth[0]/2., gx[:-1] + binwidth/2., gx[-1] + binwidth[-1]/2))

    # NOT DONE: Shift bins so the interval is '( ]' instead of '[ )'.
 
    # Determine which bin each y value belongs to. 
    binind = np.digitize(x, gx)
    # Take the sum of each bin and then divide it by the total count in each bin (average). 
    bin_means = [y[binind == i].mean() for i in range(1, len(gx))]
    bin_std = [y[binind == i].std() for i in range(1, len(gx))]
    N = [len(y[binind == i]) for i in range(1, len(gx))]

    return bin_means, bin_std, N

def bin_n(x, y, n):
    """ bin_n bins the input data 'x' into 'n' bins and calculates the mean of 'y' in those bins.
    """
    # Make sure the format is correct.
    x = standardarray(x)
    y = standardarray(y)

    # Sort the x values given.
    ind = np.argsort(x)
    x = x[ind]
    # Sort the y values according to this x-sorting.
    y = y[ind]

    # Create the bins. 
    min_i = np.arange(0,len(x),n)
    max_i = np.arange(0,len(x),n) - 1
    min_i = min_i[:-1]
    max_i = max_i[1:]
    max_i[-1] = len(x)

    # Initialize the vectors holding the binned values. 
    X = np.ones(min_i.shape)
    Y = np.ones(min_i.shape)
    S = np.ones(min_i.shape)
    Sm = np.ones(min_i.shape)
    
    # Calculate the average within each bin and put into the initialized vectors. 
    for i in range(len(max_i)):
        X[i] = np.mean(x[min_i[i]:max_i[i]])
        Y[i] = np.mean(y[min_i[i]:max_i[i]])
        S[i] = np.std(y[min_i[i]:max_i[i]])
        Sm[i] = np.std(y[min_i[i]:max_i[i]])/np.sqrt(len(y[min_i[i]:max_i[i]]))

    return X, Y, S, Sm
