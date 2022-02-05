import dadi
import numpy as np


def model_func(params, ns, pts):
    """
    Demographic model with asymmetric migrations for two populations of
    butterflies. Data and model are from McCoy et al. 2013.
    Ancestral population split into two new formed populations
    with following continuous migrations between them.

    :param nuW: Size of first new formed population.
    :param nuC: Size of second new formed population.
    :param T: Time of ancestral population split.
    :param m12: Migration rate from second population to first one.
    :param m21: Migration rate from first population to second one.
    """
    nuW, nuC, T, m12, m21 = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nuW, nu2=nuC, m12=m12, m21=m21)

    fs = dadi.Spectrum.from_phi(phi, ns, [xx, xx])
    return fs
	
	
