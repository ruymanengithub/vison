
""" 

Support Mathematical functions.

Created on Thu Aug 27 17:01:00 2020

:author: Ruyman Azzollini

""" 


# IMPORT STUFF
import numpy as np
# END IMPORT


SigmaclipResult = namedtuple('SigmaclipResult', ('clipped', 'lower', 'upper'))

def sigmaclip(a, low=4., high=4.):
    """
    Perform iterative sigma-clipping of array elements.
    
    Adapted from scipy.stats.sigmaclip. Now handling arrays with nans.

    Starting from the full sample, all elements outside the critical range are
    removed, i.e. all elements of the input array `c` that satisfy either of
    the following conditions::

        c < mean(c) - std(c)*low
        c > mean(c) + std(c)*high

    The iteration continues with the updated sample until no
    elements are outside the (updated) range.

    Parameters
    ----------
    a : array_like
        Data array, will be raveled if not 1-D.
    low : float, optional
        Lower bound factor of sigma clipping. Default is 4.
    high : float, optional
        Upper bound factor of sigma clipping. Default is 4.

    Returns
    -------
    clipped : ndarray
        Input array with clipped elements removed.
    lower : float
        Lower threshold value use for clipping.
    upper : float
        Upper threshold value use for clipping.

    Examples
    --------
    >>> from scipy.stats import sigmaclip
    >>> a = np.concatenate((np.linspace(9.5, 10.5, 31),
    ...                     np.linspace(0, 20, 5)))
    >>> fact = 1.5
    >>> c, low, upp = sigmaclip(a, fact, fact)
    >>> c
    array([  9.96666667,  10.        ,  10.03333333,  10.        ])
    >>> c.var(), c.std()
    (0.00055555555555555165, 0.023570226039551501)
    >>> low, c.mean() - fact*c.std(), c.min()
    (9.9646446609406727, 9.9646446609406727, 9.9666666666666668)
    >>> upp, c.mean() + fact*c.std(), c.max()
    (10.035355339059327, 10.035355339059327, 10.033333333333333)

    >>> a = np.concatenate((np.linspace(9.5, 10.5, 11),
    ...                     np.linspace(-100, -50, 3)))
    >>> c, low, upp = sigmaclip(a, 1.8, 1.8)
    >>> (c == np.linspace(9.5, 10.5, 11)).all()
    True

    """
    c = np.asarray(a).ravel()
    delta = 1
    while delta:
        c_std = np.nanstd(c)
        c_mean = np.nanmean(c)
        size = c.size
        critlower = c_mean - c_std * low
        critupper = c_mean + c_std * high
        c = c[(c >= critlower) & (c <= critupper)]
        delta = size - c.size

    return SigmaclipResult(c, critlower, critupper)
