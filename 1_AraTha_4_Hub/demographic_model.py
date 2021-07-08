import moments

import numpy as np

def model_func(params, ns):
    """
    Three epoch model from Huber et al., 2018.
    First epoch is ancestral.

    :param N1: Size of population during second epoch.
    :param T1: Time of second epoch.
    :param N2: Size of population during third epoch.
    :param T2: Time of third epoch.
    """
    N1, T1, N2, T2 = params

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0])
    fs = moments.Spectrum(sts)
    fs.integrate([N1], T1)
    fs.integrate([N2], T2)
    return fs
