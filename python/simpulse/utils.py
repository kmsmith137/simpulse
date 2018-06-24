import numpy as np


def make_stairsteps(data, t0, t1):
    """
    make_stairsteps(data, t0, t1) -> (tvec, yvec)

    A utility function for plotting a discrete time series as "stairsteps", in order
    to visually represent the discrete sampling.

    The 'data' argument is a 1D time series.  The 't0' and 't1' arguments are the endpoints of
    the time interval being sampled.  (More precisely, t0 is the _starting_ time of the first
    sample in the time series, and t1 is the _ending_ time of the last sample in the time series.)

    Returns a 2-tuple (tvec, yvec), suitable for plotting with matplotlib plot().

    Example::
    
       (tvec, yvec) = make_stairsteps(data, t0, t1)
       matplotlib.pyplot.plot(tvec, yvec, color='blue')
    """

    data = np.array(data, copy=False)    # convert e.g. list -> numpy array

    assert t0 < t1
    assert data.ndim == 1
    
    nt = len(data)
    tvec = np.zeros(2*nt)
    yvec = np.zeros(2*nt)

    t_delim = np.linspace(t0, t1, nt+1)
    tvec[::2] = t_delim[:-1]
    tvec[1::2] = t_delim[1:]
    yvec[::2] = data
    yvec[1::2] = data

    return (tvec, yvec)
