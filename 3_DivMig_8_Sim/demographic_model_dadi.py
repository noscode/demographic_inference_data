import dadi
import numpy


def model_func(params, ns, pts):
    """
    Three populations demographic history with small number of parameters.
    In the model ancestral population is split into population 1 and
    population 2, each of which had constant size till now days. Population 3
    is formed by split from population 2 without change of its size and had
    constant size till now too. Migration rates are symmetrical.

    :param nu1: Size of population 1.
    :param nu2: Size of population 2.
    :param nu3: Size of population 3 after split from population 2.
    :param m12: Migration rate between population 1 and population 2.
    :param m13: Migration rate between population 1 and population 3.
    :param m23: Migration rate between population 2 and population 3.
    :param T1: Time between ancestral population split and divergence of
               population 3 from population 2.
    :param T2: Time of population 3 divergence from population 2.
    """
    nu1, nu2, nu3, m12, m13, m23, T1, T2 = params
    m112 = 0

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    phi = dadi.Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nu2, m12=m112, m21=m112)

    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
    phi = dadi.Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m12, m21=m12, m23=m23, m32=m23, m13=m13, m31=m13)

    fs = dadi.Spectrum.from_phi(phi, ns, [xx, xx, xx])
    return fs
