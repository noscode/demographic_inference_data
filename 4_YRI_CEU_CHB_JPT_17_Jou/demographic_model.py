import moments
import numpy as np

def model_func(params, ns):
    """
    Demographic model for four modern human populations: YRI, CEU, CHB and JPT.
    Model and data from Jouganous et al. 2019.

    Model with sudden growth of ancestral population size, followed by split
    into population YRI and common population of CEU, CHB and JPT,
    which experience bottleneck and split with exponential recovery of all
    populations - first formation of CEU population followed by split of
    ancestal population of CHB and JPT.
    Migrations between populations are symmetrical.

    :param nuAf: The ancestral population size after sudden growth
                 and size of YRI population.
    :param nuB: The bottleneck size of CEU+CHB common population.
    :param nuEu0: The bottleneck size for CEU population.
    :param nuEu: The final size of CEU population after exponential growth.
    :param nuAs0: The bottleneck size for CHB population.
    :param nuAs: The final size of CHB population after exponential growth.
    :param nuJp0: The bottleneck size for JPT population.
    :param nuJp: The final size of JPT population after exponential growth.
    :param mAfB: The scaled symmetric migration rate between YRI and CEU+CHB  
                 populations.
    :param mAfEu: The scaled symmetric migration rate between YRI and CEU
                  populations.
    :param mAfAs: The scaled symmetric migration rate between YRI and CHB
                  populations.
    :param mEuAs: The scaled symmetric migration rate between CEU and CHB
                  populations.
    :param mChJp: The scaled symmetric migration rate between CHB and JPT
                  populations.
    :param TAf: The scaled time between ancestral population growth
                and first split.
    :param TB: The time between the first split and second. Time of CEU+CHB+JPT
               population existence.
    :param TEuAs: The time between second split and third split. Time of
                  CHB+JPT population existence.
    :param TEuAs: The time between third split and present.
    """
    nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs, nuJp0, nuJp, mAfB, mAfEu, mAfAs, mEuAs, mChJp, TAf, TB, TEuAs, TChJp = params

    sts = moments.LinearSystem_1D.steady_state_1D(sum(ns))
    fs = moments.Spectrum(sts)

    fs.integrate([nuAf], TAf)
    
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], sum(ns[1:]))
    
    m = np.array([[0, mAfB],[mAfB, 0]])
    fs.integrate([nuAf, nuB], TB, m=m)
    
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], sum(ns[2:]))

    nuEu_func = lambda t: nuEu0 * (nuEu / nuEu0) ** (t / (TEuAs + TChJp))
    nuAs_func = lambda t: nuAs0 * (nuAs / nuAs0) ** (t / (TEuAs + TChJp))
    nu_func = lambda t: [nuAf, nuEu_func(t), nuAs_func(t)]
    m = np.array([[0, mAfEu, mAfAs],
                  [mAfEu, 0, mEuAs],
                  [mAfAs, mEuAs, 0]])
    fs.integrate(nu_func, TEuAs, m=m)

    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], ns[3])

    nuEu_func = lambda t: nuEu0 * (nuEu / nuEu0) ** ((t + TEuAs) / (TEuAs + TChJp))
    nuAs_func = lambda t: nuAs0 * (nuAs / nuAs0) ** ((t + TEuAs) / (TEuAs + TChJp))
    nuJp_func = lambda t: nuJp0 * (nuJp / nuJp0) ** (t / (TChJp))
    nu_func = lambda t: [nuAf, nuEu_func(t), nuAs_func(t), nuJp_func(t)]
    m = np.array([[0, mAfEu, mAfAs, 0],
                  [mAfEu, 0, mEuAs, 0],
                  [mAfAs, mEuAs, 0, mChJp],
                  [0, 0, mChJp, 0]])
    fs.integrate(nu_func, TChJp, m=m)

    return fs

