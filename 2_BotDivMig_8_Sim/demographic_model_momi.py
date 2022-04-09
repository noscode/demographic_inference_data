import momi
import numpy as np

def model_func(params):
    Nanc, nu_gen, f, nu1_gen, nu2_gen, m12_gen, m21_gen, T1_gen, T2_gen = params

    nu1_gen_init = nu_gen * f
    nu2_gen_init = nu_gen * (1 - f)

    model = momi.DemographicModel(N_e=1)
    model.add_leaf("pop0", N=nu1_gen*Nanc, g=np.log(nu1_gen/nu1_gen_init) / (T2_gen*2*Nanc))
    model.add_leaf("pop1", N=nu2_gen*Nanc, g=np.log(nu2_gen/nu2_gen_init) / (T2_gen*2*Nanc))

    model.move_lineages("pop1", "pop0", t=T2_gen*2*Nanc, N=nu_gen*Nanc, g=0)
    model.set_size("pop0", N=Nanc, g=0, t=(T1_gen+T2_gen)*2*Nanc)
    return model
