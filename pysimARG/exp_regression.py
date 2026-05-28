import numpy as np
from scipy.optimize import curve_fit


def shifted_exponential(x, a, b, c):
    """
    a: Amplitude
    b: Growth/decay rate
    c: Vertical shift (asymptote)
    """
    return a * np.exp(-b * x) + c


def exp_regression(x, y):
    """
    Fit an exponential regression model to the given data.
    
    The model is defined as:
        y = a * exp(-b * x) + c
    
    Parameters
    ----------
    x : numpy array
        Independent variable data.
    y : numpy array
        Dependent variable data.
    Returns
    -------
    tuple
        Optimal parameters (a, b, c) for the fitted model.
    """
    initial_guess = [1.0, 1.0, np.min(y)] 
    popt, pcov = curve_fit(shifted_exponential, x, y, p0=initial_guess)
    return tuple(popt)
