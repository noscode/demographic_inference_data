import moments
import numpy as np


def model_func(params, ns):
    """
    Four population demographic history with 18 parameters. Ancestral
    population of constant size was split (T1 + T2 + T3) time ago into two new
    populations. First (1) population has constant size till nowdays.
    The second formed population turned out to be a common population to
    another three populations: (T2 + T3) time ago it splits and formed
    so-called population 2 and T3 time ago the second formed population
    divided into population 3 and population 4. All migrations between
    populations were symmetrical.    

    :param nu1: Size of population 1 after split of ancestral population.
    :param nu234: Size of common ancestor population of populations 2,
                  3 and 4 after split of ancestral population.
    :param nu2: Size of population 2 after split of common ancestor population
                of populations 2, 3 and 4.
    :param nu34: Size of common ancestor population of populations 3 and 4
                 after division of population 2 from their common ancestor
                 population.
    :param nu3: Size of population 3.
    :param nu4: Size of population 4.
    :param m1_234: Symmetric migration rate between population 1 and common
                   ancestor population of 2, 3 and 4 populations (between first
                   and second splits).
    :param m1_2: Symmetric migration rate between population 1 and population
                 2 (between second and third splits).
    :param m1_34: Symmetric migration rate between population 1 and common
                  ancestor of population 3 and population 4 (between second and
                  third splits).
    :param m2_34: Symmetric migration rate between population 2 and common
                  ancestor of population 3 and population 4 (between second and
                  third splits).
    :param m1_3: Symmetric migration rate between population 1 and population 3
                 (after the third split).
    :param m1_4: Symmetric migration rate between population 1 and population 4
                 (after the third split).
    :param m2_3: Symmetric migration rate between population 2 and population 3
                 (after the third split).
    :param m2_4: Symmetric migration rate between population 2 and population 4
                 (after the third split).
    :param m3_4: Symmetric migration rate between population 3 and population 4
                 (after the third split).
    :param T1: Time between ancestral population split (population 1 formation)
               and next split.
    :param T2: Time between ancestral population of populations 2, 3 and 4
               split (population 2 formation) and next split.
    :param T3: Time of ancestral population of populations 3 and 4 split which
               have led to formations of population 3 and population 4.
    """
    nu1, nu234, nu2, nu34, nu3, nu4, m1_234, m1_2, m1_34, m2_34, m1_3, m1_4, m2_3, m2_4, m3_4, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(sum(ns))
    fs = moments.Spectrum(sts)

    fs = moments.Manips.split_1D_to_2D(fs, ns[0], sum(ns[1:]))
    m = np.array([
            [0,     m1_234], 
            [m1_234, 0]
            ])
    fs.integrate(Npop=[nu1, nu234], tf=T1, m=m)
    
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], sum(ns[2:]))
    m = np.array([
            [0,     m1_2,     m1_34], 
            [m1_2,     0,     m2_34],
            [m1_34, m2_34,     0]
            ])
    fs.integrate(Npop=[nu1, nu2, nu34], tf=T2, m=m)

    fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], sum(ns[3:]))
    m = np.array([
            [0,     m1_2,     m1_3,     m1_4], 
            [m1_2,     0,     m2_3,     m2_4],
            [m1_3,     m2_3,     0,     m3_4],
            [m1_4,     m2_4,     m3_4,    0]
            ])
    fs.integrate(Npop=[nu1, nu2, nu3, nu4], tf=T3, m=m)

    return fs
    
    
