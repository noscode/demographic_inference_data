import moments
import numpy as np

def model_func(params, ns):
    """
    Demographic history of five populations from Ragsdale and Gravel 2019.

    This model is part of the stdpopsim models for Homo Sapiens.
    See more information (OutOfAfricaArchaicAdmixture_5R19 model):
    https://popsim-consortium.github.io/stdpopsim-docs/stable/catalog.html

    The code is based on the examples of dadi for CUDA:
    https://bitbucket.org/gutenkunstlab/dadi/src/master/examples/CUDA/models_moments.py

    The model describes the splits of neandertalian and archaic african
    populations from ancestral population. Then the history of this population
    corresponds to the usual Out-of-Africa three population history.
    At the end two archaic populations are removed and final spectrum
    is three dimentional.

    :param nu_YRI: The African population size after sudden growth
                   and size of YRI population.
    :param nu_B: The bottleneck size of CEU+CHB common population.
    :param nu_CEU0: The bottleneck size for CEU population.
    :param nu_CEU: The final size of CEU population after exponential growth.
    :param nu_CHB0: The bottleneck size for CHB population.
    :param nu_CHB: The final size of CHB population after exponential growth.
    :param m_AF_B: The scaled symmetric migration rate between YRI and CEU+CHB  
                 populations.
    :param m_YRI_CEU: The scaled symmetric migration rate between YRI and CEU
                  populations.
    :param m_YRI_CHB: The scaled symmetric migration rate between YRI and CHB
                  populations.
    :param m_CEU_CHB: The scaled symmetric migration rate between CEU and CHB
                  populations.
    :param m_AF_arch_af: The scaled symmetric migration rate between YRI and
                  archaic African populations that begins after sudden growth
                  of YRI population.
    :param m_OOA_nean: The scaled symmetric migration rate between
                  Out of AFrica populations and
                  neanderthalian population.
    :param T_Nean_arch_AF: The time between neanderthalian split and
               archaic african population split.
    :param T_AF_arch_AF: The time between archaic african population split
               and sudden growth of YRI population.
    :param T_AF_no_arch_mig: The time between YRI sudden population growth
               and start of admixture from archaic african population.
    :param T_AF_arch_mig: The time between start of admixture from archaic
               african population and OOA event.
    :param T_Bot: The time between the OOA event and Eurasian population split.
               Time of CEU+CHB population existence.
    :param T_EU_AS_arch_adm_end: The time between Eurasian population split and
               before end of archaic admixture.
    :param T_EU_AS_no_arch_adm: The time between end of archaic admixture and
               present.
    """
    nu_YRI, nu_B, nu_CEU0, nu_CEU, nu_CHB0, nu_CHB, m_AF_B, m_YRI_CEU, m_YRI_CHB, m_CEU_CHB, m_AF_arch_af, m_OOA_nean, T_Nean_arch_AF,  T_AF_arch_AF, T_AF_no_arch_mig,  T_AF_arch_mig, T_Bot, T_EU_AS_arch_adm_end, T_EU_AS_no_arch_adm = params

    n1, n2, n3, n4, n5 = ns

    # We want an order of Nean, Archaic_Afr, YRI, CEU, CHB

    sts = moments.LinearSystem_1D.steady_state_1D(n1+n2+n3+n4+n5)
    fs = moments.Spectrum(sts)

    # Split off Neanderthal pop
    # [Nean, others]
    fs = moments.Manips.split_1D_to_2D(fs, n1, n2+n3+n4+n5)
    fs.integrate([1., 1.],  T_Nean_arch_AF)
    # Split off archaic African pop
    # [Nean, archaic, others]
    fs = moments.Manips.split_2D_to_3D_2(fs, n2, n3+n4+n5)
    fs.integrate([1., 1., 1.],  T_AF_arch_AF)
    # African population growth
    # [Nean, archaic, YRI]
    fs.integrate([1., 1., nu_YRI],  T_AF_no_arch_mig)
    # Archaic African migration begins
    # migratiion between archaic and YRI
    mig_mat = [[0, 0, 0],
               [0, 0, m_AF_arch_af],
               [0, m_AF_arch_af, 0]]
    fs.integrate([1., 1., nu_YRI], T_AF_arch_mig, m=mig_mat)
    # Split of Eurasian ancestral pop
    # [Nean, archaic, YRI, others]
    # migrations berween archaic and YRI, Nean and others, Af and B
    fs = moments.Manips.split_3D_to_4D_3(fs, n3, n4+n5)
    nus = [1.0, 1.0, nu_YRI, nu_B]
    mig_mat = [[0, 0, 0, m_OOA_nean],
               [0, 0, m_AF_arch_af, 0],
               [0, m_AF_arch_af, 0, m_AF_B],
               [m_OOA_nean, 0, m_AF_B, 0]]
    fs.integrate(nus, T_Bot, m=mig_mat)
    # Split of European and Asian ancestral pops
    # [Nean, archaic, YRI, CEU, CHB]
    fs = moments.Manips.split_4D_to_5D_4(fs, n4, n5) 
    nuEu_func = lambda t: nu_CEU0*(nu_CEU/nu_CEU0)**(t/(T_EU_AS_arch_adm_end + T_EU_AS_no_arch_adm))
    nuAs_func = lambda t: nu_CHB0*(nu_CHB/nu_CHB0)**(t/(T_EU_AS_arch_adm_end + T_EU_AS_no_arch_adm))
    nu_func = lambda t: [1.0, 1.0, nu_YRI, nuEu_func(t), nuAs_func(t)]
    mig_mat = [[0, 0, 0, m_OOA_nean, m_OOA_nean],
               [0, 0, m_AF_arch_af, 0, 0],
               [0, m_AF_arch_af, 0, m_YRI_CEU, m_YRI_CHB],
               [m_OOA_nean, 0, m_YRI_CEU, 0, m_CEU_CHB],
               [m_OOA_nean, 0, m_YRI_CHB, m_CEU_CHB, 0]]
    fs.integrate(nu_func, T_EU_AS_arch_adm_end, m=mig_mat)
    # End of archaic migration. Remove archaic pops for efficiency.
    fs = fs.marginalize([0, 1])

    _, _, _, nuEu_temp, nuAs_temp = nu_func(T_EU_AS_arch_adm_end)
    nu_func = lambda t: [nu_YRI,
                         nuEu_temp * (nu_CEU/nuEu_temp)**(t/T_EU_AS_no_arch_adm),
                         nuAs_temp * (nu_CHB/nuAs_temp)**(t/T_EU_AS_no_arch_adm)]

    # We use initial_t argument so we can reuse nuEu_func and nuAs_func
    mig_mat = [[0, m_YRI_CEU, m_YRI_CHB],
               [m_YRI_CEU, 0, m_CEU_CHB],
               [m_YRI_CHB, m_CEU_CHB, 0]]
    fs.integrate(nu_func, T_EU_AS_no_arch_adm, m=mig_mat)

    return fs

