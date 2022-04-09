import momi
import numpy as np

def model_func(params):
    Nanc, nuB_gen, nuF_gen, tB_gen, tF_gen = params

    model = momi.DemographicModel(N_e=1)
    model.add_leaf("pop0", N=nuF_gen*Nanc, g=0)
    model.set_size("pop0", N=nuB_gen*Nanc, g=0, t=tF_gen*2*Nanc)
    model.set_size("pop0", N=Nanc, g=0, t=(tF_gen+tB_gen)*2*Nanc)
    return model
