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

    sts = moments.LinearSystem_1D.steady_state_1D(ns[0])
    fs = moments.Spectrum(sts)
    nu1_func = lambda t: [1.0 * (nu1 / 1.0) ** (t / t2)]
    fs.integrate(nu1_func, t1)
    nu1_func = lambda t: [nu1 * (nu2 / nu1) ** (t / t2)]
    fs.integrate(nu1_func, t2)
    nu1_func = lambda t: [nu2 * (nu3 / nu2) ** (t / t3)]
    fs.integrate(nu1_func, t3)
    nu1_func = lambda t: [nu3 * (nu4 / nu3) ** (t / t4)]
    fs.integrate(nu1_func, t4)
    nu1_func = lambda t: [nu4 * (nu5 / nu4) ** (t / t5)]
    fs.integrate(nu1_func, t5)
    fs.integrate([nu5], t6)
    return fs
