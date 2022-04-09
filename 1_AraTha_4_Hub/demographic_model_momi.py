import momi
import numpy as np

def model_func(params):
    Nanc, N1_gen, T1_gen, N2_gen, T2_gen = params

    model = momi.DemographicModel(N_e=1)
    model.add_leaf("pop0", N=N2_gen*Nanc, g=0)
    model.set_size("pop0", N=N1_gen*Nanc, g=0, t=T1_gen*2*Nanc)
    model.set_size("pop0", N=Nanc, g=0, t=(T1_gen+T2_gen)*2*Nanc)
    return model
