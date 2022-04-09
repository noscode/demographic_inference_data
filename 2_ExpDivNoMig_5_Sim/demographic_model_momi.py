import momi
import numpy as np

def model_func(params):
    Nanc, nu_gen, nu1_gen, nu2_gen, T1_gen, T2_gen = params

    model = momi.DemographicModel(N_e=1)
    model.add_leaf("pop0", N=nu1_gen*Nanc, g=0)
    model.add_leaf("pop1", N=nu2_gen*Nanc, g=0)

    model.move_lineages("pop0", "pop1", t=T2_gen*2*Nanc, N=Nanc, g=np.log(nu_gen/1.0) / (T1_gen*2*Nanc))
    nodel.set_size("pop0", N=nu_gen*Nanc, g=0, t=(T1_gen+T2_gen)*2*Nanc)
    return model
