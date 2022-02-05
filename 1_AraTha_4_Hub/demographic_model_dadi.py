import dadi

import numpy as np

def model_func(params, ns, pts):
    """
    Three epoch model from Huber et al., 2018.
    First epoch is ancestral.

    :param N1: Size of population during second epoch.
    :param T1: Time of second epoch.
    :param N2: Size of population during third epoch.
    :param T2: Time of third epoch.
    """
    N1, T1, N2, T2 = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi, xx, T1, nu=N1)
    phi = dadi.Integration.one_pop(phi, xx, T2, nu=N2)

    fs = dadi.Spectrum.from_phi(phi, ns, [xx])
    return fs
