import dadi
import numpy as np


def model_func(params, ns, pts):
    """
    Demographic model of isolation for two populations with exponential growth
    of an ancestral population followed by split.

    :param nu: Size of ancestral population after exponential growth.
    :param nu1: Size of population 1 after split.
    :param nu2: Size of population 2 after split.
    :param T1: Time between exponential growth of ancestral population and its
               split.
    :param T2: Time of ancestral population split.
    """
    nu, nu1, nu2, T1, T2 = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    nu_func = lambda t: (nu) ** (t / T1)
    phi = dadi.Integration.one_pop(phi, xx, T1, nu=nu_func)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=nu1, nu2=nu2)

    fs = dadi.Spectrum.from_phi(phi, ns, [xx, xx])
    return fs
    
    
