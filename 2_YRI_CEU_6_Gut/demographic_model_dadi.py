import dadi
import numpy


def model_func(params, ns, pts):
    """
    Demographic model for two modern human populations: YRI and CEU.
    Data and model are from Gutenkunst et al., 2009.

    Model with sudden growth of ancestral population size, followed by split,
    bottleneck in second population (CEU) with exponential recovery and
    symmetric migration.

    :param nu1F: The ancestral population size after growth.
    :param nu2B: The bottleneck size for second population (CEU).
    :param nu2F: The final size for second population (CEU).
    :param m: The scaled symmetric migration rate.
    :param Tp: The scaled time between ancestral population growth
               and the split.
    :param T: The time between the split and present.
    """
    nu1F, nu2B, nu2F, m, Tp, T = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi, xx, Tp, nu=nu1F)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    nu2_func = lambda t: nu2B*(nu2F/nu2B)**(t/T)
    phi = dadi.Integration.two_pops(phi, xx, T, nu1=nu1F, nu2=nu2_func, m12=m, m21=m)

    fs = dadi.Spectrum.from_phi(phi, ns, [xx, xx])
    return fs
	
	
