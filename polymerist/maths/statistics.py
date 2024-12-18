'''For computing statistical measures and series' of structured data'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import numpy as np
from dataclasses import dataclass


# CALCULATION OF STATISTICS
def RMSE(pred : np.ndarray, obs : np.ndarray) -> float:
    '''Computes root-mean squared error between predicted and observed sets of values'''
    sq_err = (obs - pred)**2
    return np.sqrt(sq_err.mean())

@dataclass
class Accumulator:
    '''Compact container for accumulating averages'''
    sum : float = 0.0
    count : int = 0

    @property
    def average(self) -> float:
        return self.sum / self.count
    
    
# SERIES-TO-SERIES OPERATIONS
def normalize(series : np.ndarray) -> np.ndarray:
    '''Normalize a series of data to have all values lie in the range [-1, 1]'''
    return (series - series.min()) / (series.max() - series.min())

def standardize(series : np.ndarray) -> np.ndarray:
    '''Standardize a series of data to have mean 0 and standard deviation of 1'''
    return (series - series.mean()) / series.std() 

def autocorrelate(series : np.ndarray) -> np.ndarray:
    '''Compute autocorrelation of a vector series of data points'''
    assert(series.ndim == 1) # only operate for vectors of data
    
    series = standardize(series)
    autocorr = np.correlate(series, series, mode='full') # need full mode to get proper length
    autocorr = autocorr[autocorr.size//2:] # only keep second half of array (Hermitian, so symmetric about midpoint)
    autocorr /= series.size # normalize to range [-1, 1] by including number of datapoints (due to variance estimator)

    return autocorr
