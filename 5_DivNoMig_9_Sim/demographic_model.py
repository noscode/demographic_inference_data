import moments
import numpy as np


def model_func(params, ns):
    """
    Simple demographic history of five populations without migrations.
    Model have 9 parameters. Each population number i has constant size
    of nui. Ancestral population splits (T1+T2+T3+T4) time ago to
    population 1 and 2. Then (T2+T3+T4) time ago population 2 splits into
    population 2 and 3. (T3 + T4) time ago population 3 split to
    populations 3 and 4. And finally T4 time ago population 4 split
    in populations 4 and 5.

    :param nu1: Size of population 1.
    :param nu2: Size of population 2.
    :param nu3: Size of population 3 after split from population 2.
    :param nu4: Size of population 4 after split from population 3.
    :param nu5: Size of population 5 after split from population 4.
    :param T1: Time between ancestral population split and second split.
    :param T2: Time between second and third splits.
    :param T3: Time between third and fourth splits.
    :param T4: Time of fourth split.
    """
    nu1, nu2, nu3, nu4, nu5, T1, T2, T3, T4 = params
    sts = moments.LinearSystem_1D.steady_state_1D(sum(ns))
    fs = moments.Spectrum(sts)

    fs = moments.Manips.split_1D_to_2D(fs, ns[0], sum(ns[1:]))

    fs.integrate(Npop=[nu1, nu2], tf=T1)
    
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], sum(ns[2:]))
    
    fs.integrate(Npop=[nu1, nu2, nu3], tf=T2)

    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], sum(ns[3:]))
    
    fs.integrate(Npop=[nu1, nu2, nu3, nu4], tf=T3)

    fs = moments.Manips.split_4D_to_5D_4(fs, ns[3], ns[4])

    fs.integrate(Npop=[nu1, nu2, nu3, nu4, nu5], tf=T4)

    return fs
    
    
