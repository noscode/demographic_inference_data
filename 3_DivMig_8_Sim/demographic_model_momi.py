import momi
import numpy as np

def model_func(params):
    Nanc, nu1_gen, nu2_gen, nu3_gen, T1_gen, T2_gen = params

    model = momi.DemographicModel(N_e=1)
    model.add_leaf("pop0", N=nu1_gen*Nanc, g=0)
    model.add_leaf("pop1", N=nu2_gen*Nanc, g=0)
    model.add_leaf("pop2", N=nu3_gen*Nanc, g=0)

    model.move_lineages("pop2", "pop1", t=T2_gen*2*Nanc, N=nu2_gen*Nanc, g=0)
    model.move_lineages("pop1", "pop0", t=(T1_gen+T2_gen)*2*Nanc, N=Nanc, g=0)
    return model
