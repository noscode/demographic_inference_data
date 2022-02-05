import moments

import numpy as np

def model_func(params, ns):
    """
    ZigZag model from Stephan and Durbin, 2014.

    :param nu1: Size of population after first exponential growth.
    :param nu2: Size of population after first exponential decrease.
    :param nu3: Size of population after second exponential growth.
    :param nu4: Size of population after second exponential decrease.
    :param nu5: Size of population now after third exp. growth.
    :param t1: Time of first exponential growth.
    :param t2: Time of first exponential decrease.
    :param t3: Time of second exponential growth.
    :param t4: Time of second exponential decrease.
    :param t5: Time of third exponential growth.
    :param t6: Time of end of the last growth.
    """
    nu1, nu2, nu3, nu4, nu5, t1, t2, t3, t4, t5, t6 = params


    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    nu1_func = lambda t: 1.0 * (nu1 / 1.0) ** (t / t1)
    phi = dadi.Integration.one_pop(phi, xx, t1, nu=nu1_func)
    nu1_func = lambda t: nu1 * (nu2 / nu1) ** (t / t2)
    phi = dadi.Integration.one_pop(phi, xx, t2, nu=nu1_func)
    nu1_func = lambda t: nu2 * (nu3 / nu2) ** (t / t3)
    phi = dadi.Integration.one_pop(phi, xx, t3, nu=nu1_func)
    nu1_func = lambda t: nu3 * (nu4 / nu3) ** (t / t4)
    phi = dadi.Integration.one_pop(phi, xx, t4, nu=nu1_func)
    nu1_func = lambda t: nu4 * (nu5 / nu4) ** (t / t5)
    phi = dadi.Integration.one_pop(phi, xx, t5, nu=nu1_func)
    phi = dadi.Integration.one_pop(phi, xx, t6, nu=nu5)

    fs = dadi.Spectrum.from_phi(phi, ns, [xx])
    return fs
