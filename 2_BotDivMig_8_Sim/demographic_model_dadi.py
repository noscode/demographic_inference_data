import dadi
import numpy as np


def model_func(params, ns, pts):
    """
    Demographic history of two populations with bottleneck of ancestral
    population followed by split and growth of both new formed populations
    exponentially and linearly correspondingly.

    :param nu: Size of ancestral population after sudden decrease.
    :param f: Fraction in which ancestral population splits.
    :param nu1: Size of population 1 after exponential growth.
    :param nu2: Size of population 2 after linear growth.
    :param m12: Migration rate from subpopulation 2 to subpopulation 1.
    :param m21: Migration rate from subpopulation 1 to subpopulation 2.
    :param T1: Time between sudden growth of ancestral population and its
               split.
    :param T2: Time of ancestral population split.
    """
    nu, f, nu1, nu2, m12, m21, T1, T2 = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi, xx, T1, nu=nu)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nu1_init = nu * f
    nu2_init = nu * (1 - f)
    nu1_func = lambda t: nu1_init + (nu1 - nu1_init) * (t / T2)
    nu2_func = lambda t: nu2_init * (nu2 / nu2_init) ** (t / T2)
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=nu1_func, nu2=nu2_func, m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, [xx, xx])
    return fs
    
    
