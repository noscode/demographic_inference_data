import momi
import numpy as np

def model_func(params):
    Nanc, N1_gen, T1_gen, N2_gen, T2_gen = params

    model = momi.DemographicModel(N_e=1e5)
    model.add_leaf("Pop1", N=N2_gen*Nanc, g=0)
    model.set_size("Pop1", N=N1_gen*Nanc, g=0, t=T1)
    model.set_size("Pop1", N=Nanc, g=0, t=T1+T2)
    return model
