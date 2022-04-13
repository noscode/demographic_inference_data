import momi
import numpy as np

def model_func(params):
    Nanc, nu1F_gen, nu2B_gen, nu2F_gen, Tp_gen, T_gen = params

    model = momi.DemographicModel(N_e=1)
    model.add_leaf("pop0", N=nu1F_gen*Nanc, g=0)
    model.add_leaf("pop1", N=nu2F_gen*Nanc, g=np.log(nu2F_gen/nu2B_gen) / (T_gen*2*Nanc))

    model.move_lineages("pop1", "pop0", t=T_gen*2*Nanc, N=nu1F_gen*Nanc, g=0)
    model.set_size("pop0", N=Nanc, g=0, t=(Tp_gen+T_gen)*2*Nanc)

    return model
