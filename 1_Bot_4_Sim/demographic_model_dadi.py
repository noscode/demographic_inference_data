import dadi

import numpy as np

def model_func(params, ns, pts):
    """
    Classical one population bottleneck model.

    :param nuB: Size of population during bottleneck.
    :param nuF: Size of population now.
    :param tB: Time of bottleneck duration.
    :param tF: Time after bottleneck finished.
    """
    nuB, nuF, tB, tF = params

    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)

    phi = dadi.Integration.one_pop(phi, xx, tB, nu=nuB)
    phi = dadi.Integration.one_pop(phi, xx, tF, nu=nuF)

    fs = dadi.Spectrum.from_phi(phi, ns, [xx])
    return fs
