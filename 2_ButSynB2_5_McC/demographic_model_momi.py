import momi
import numpy as np

def model_func(params):
    Nanc, nuW_gen, nuC_gen, T_gen = params

    model = momi.DemographicModel(N_e=1)
    model.add_leaf("WY", N=nuW_gen*Nanc, g=0)
    model.add_leaf("CO", N=nuC_gen*Nanc, g=0)

    model.move_lineages("CO", "WY", t=T_gen*2*Nanc, N=Nanc, g=0)
    return model
