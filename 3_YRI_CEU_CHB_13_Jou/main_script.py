import moments
from demographic_model import model_func

n_pop = 3
pop_labels = ["YRI", "CEU", "CHB"]

par_labels = ['nuAf', 'nuB', 'nuEu0', 'nuEu', 'nuAs0', 'nuAs',
              'mAfB', 'mAfEu', 'mAfAs', 'mEuAs', 'TAf', 'TB', 'TEuAs']

# Parameters from Jouganous et al. (2017) paper (Table 4):
# N_A = 11293
# N_Af = 24486
# N_B = 3034
# N_Eu0 = 2587
# r_Eu = 0.17e-2  # it is percent in table
# N_As0 = 958
# r_As = 0.30e-2  # it is percent in table
# m_Af_B = 15.6e-5
# m_Af_Eu = 1.00e-5
# m_Af_As = 0.48e-5
# m_Eu_As = 3.99e-5
# T_Af = 349e3  # kya
# T_B = 121e3  # kya
# T_Eu_As = 44e3  # kya
# We traslated them and got:
popt = [2.168245815992208, 0.2686620030107146, 0.22907996103781103,
        3.014506849317947, 0.08483131143186044, 7.9870350674616715,
        3.523416, 0.22586,
        0.10841279999999999, 0.9011814,
        0.3480947917080156, 0.11755832877858423, 0.0671761878734767]

lower_bound = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2,
               0, 0, 0, 0, 1e-15, 1e-15, 1e-15]
upper_bound = [100, 100, 100, 100, 100, 100, 10, 10, 10, 10, 5, 5, 5]

mu = 1.44e-8  # mutation rate
L = None  # effective length of sequence

ns = [80, 80, 80]

# Get maximum log-likelihood
model = model_func(popt, ns)
data = moments.Spectrum.from_file("fs_data.fs")
max_ll = moments.Inference.ll_multinom(model, data)

# Get ancestral population size
theta =  moments.Inference.optimal_sfs_scaling(model, data) # mutation flux
Nanc = 11293

if __name__ == "__main__":
    print('Maximum log composite likelihood: {0}'.format(max_ll))

    print('Optimal value of theta: {0}'.format(theta))

    print('Size of the ancestral population: {0}'.format(Nanc))

    # Draw model
    model = moments.ModelPlot.generate_model(model_func, popt,
                                             [4 for _ in range(n_pop)])
    moments.ModelPlot.plot_model(model,
        save_file='model_plot.png',
        fig_title=f'Demographic model for 3_YRI_CEU_CHB_13_Jou',
        pop_labels=pop_labels,
        nref=Nanc,
        draw_scale=True,
        draw_ancestors=True,
        gen_time=1.0,
        gen_time_units='Generations',
        reverse_timeline=True)
    print('Model plot is saved to model_plot.png')
