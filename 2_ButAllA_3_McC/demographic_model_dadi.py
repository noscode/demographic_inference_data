import dadi
import numpy as np


def model_func(params, ns, pts):
    """
    Demographic model without migration for two populations of butterflies.
    Data and model are from McCoy et al. 2013.
    Model is very simple: ancestral population splits into two new populations
    of constant size.

    :param nuW: Size of first subpopulation.
    :param nuC: Size of second subpopulation.
    :param T: Time of ancestral population split.
    """
    nuW, nuC, T = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nuW, nu2=nuC)

    fs = dadi.Spectrum.from_phi(phi, ns, [xx])
    return fs
	
	
