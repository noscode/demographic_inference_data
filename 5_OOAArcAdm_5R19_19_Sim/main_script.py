import moments
from demographic_model import model_func
import numpy as np

n_pop = 5
pop_labels = ['Neand', 'AA', 'YRI', 'CEU', 'CHB']

par_labels = ['nu_YRI', 'nu_B', 'nu_CEU0', 'nu_CEU', 'nu_CHB0', 'nu_CHB',
              'm_AF_B', 'm_YRI_CEU', 'm_YRI_CHB', 'm_CEU_CHB',
              'm_AF_arch_af', 'm_OOA_nean',
              'T_Nean_arch_AF',  'T_AF_arch_AF', 'T_AF_no_arch_mig',
              'T_AF_arch_mig', 'T_Bot',
              'T_EU_AS_arch_adm_end', 'T_EU_AS_no_arch_adm']

# First we set out the maximum likelihood values of the various parameters
# given in Table 1 (under archaic admixture).
N_0 = 3600
nu_YRI = 13900 / N_0
nu_B = 880 / N_0
nu_CEU0 = 2300 / N_0
nu_CHB0 = 650 / N_0

# Times are provided in years, so we convert into generations.
# In the published model, the authors used a generation time of 29 years to
# convert from genetic to physical units
generation_time = 29

T_AF = 300e3 / (generation_time * 2 * N_0)
T_B = 60.7e3 / (generation_time * 2 * N_0)
T_EU_AS = 36.0e3 / (generation_time * 2 * N_0)
T_arch_afr_split = 499e3 / (generation_time * 2 * N_0)
T_arch_afr_mig = 125e3 / (generation_time * 2 * N_0)
T_nean_split = 559e3 / (generation_time * 2 * N_0)
T_arch_adm_end = 18.7e3 / (generation_time * 2 * N_0)

# We need to work out the starting (diploid) population sizes based on
# the growth rates provided for these two populations
r_CEU = 0.00125
r_CHB = 0.00372
nu_CEU = nu_CEU0 / np.exp(-r_CEU * T_EU_AS * 2 * N_0)
nu_CHB = nu_CHB0 / np.exp(-r_CHB * T_EU_AS * 2 * N_0)

# Migration rates during the various epochs.
m_AF_B = 52.2e-5 * 2 * N_0
m_YRI_CEU = 2.48e-5 * 2 * N_0
m_YRI_CHB = 0e-5 * 2 * N_0
m_CEU_CHB = 11.3e-5 * 2 * N_0
m_AF_arch_af = 1.98e-5 * 2 * N_0
m_OOA_nean = 0.825e-5 * 2 * N_0

T_Nean_arch_AF = T_nean_split-T_arch_afr_split
T_AF_arch_AF = T_arch_afr_split-T_AF
T_AF_no_arch_mig = T_AF-T_arch_afr_mig
T_AF_arch_mig = T_arch_afr_mig-T_B
T_Bot = T_B-T_EU_AS
T_EU_AS_arch_adm_end = T_EU_AS-T_arch_adm_end
T_EU_AS_no_arch_adm = T_arch_adm_end

popt = [
    nu_YRI, nu_B, nu_CEU0, nu_CEU, nu_CHB0, nu_CHB,
    m_AF_B, m_YRI_CEU, m_YRI_CHB, m_CEU_CHB, m_AF_arch_af, m_OOA_nean,
    T_Nean_arch_AF,  T_AF_arch_AF, T_AF_no_arch_mig,  T_AF_arch_mig,
    T_Bot, T_EU_AS_arch_adm_end, T_EU_AS_no_arch_adm,
]

lower_bound = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2,
               0, 0, 0, 0, 0, 0,
               1e-15, 1e-15, 1e-15, 1e-15, 1e-15, 1e-15, 1e-15]
upper_bound = [100, 100, 100, 100, 100, 100,
               10, 10, 10, 10, 10, 10,
               5, 5, 5, 5, 5, 5, 5]

mu = 1.29e-08  # mutation rate
L = 20000000  # effective length of sequence
Nanc = N_0
theta = 4 * mu * L * Nanc  # mutation flux

ns_per_pop = 10
ns = [ns_per_pop for _ in range(n_pop)]

# Get maximum log-likelihood
model = model_func(popt, ns)
data = model * theta
max_ll = moments.Inference.ll_multinom(model, data)

if __name__ == "__main__":
    data.to_file('fs_data.fs')
    print('Simulated data saved to fs_data.fs')

    print('Maximum log composite likelihood: {0}'.format(max_ll))

    theta = moments.Inference.optimal_sfs_scaling(model, data)
    print('Optimal value of theta: {0}'.format(theta))

    theta0 = 4 * mu * L
    Nanc = int(theta / theta0)
    print('Size of the ancestral population: {0}'.format(Nanc))

    # Draw model
    model = moments.ModelPlot.generate_model(model_func, popt,
                                             [4 for _ in range(n_pop)])
    moments.ModelPlot.plot_model(model,
        save_file='model_plot.png',
        fig_title=f'5_OOAArcAdm_5R19_19_Sim, Nanc: {Nanc}',
        pop_labels=pop_labels,
        nref=Nanc,
        draw_scale=False,
        draw_ancestors=True,
        gen_time=1.0,
        gen_time_units='Generations',
        reverse_timeline=True)
    print('Model plot is saved to model_plot.png')
