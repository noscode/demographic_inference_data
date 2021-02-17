import moments
from demographic_model import model_func

n_pop = 4
pop_labels = ["YRI", "CEU", "CHB", "JPT"]
par_labels = ['nuAf', 'nuB', 'nuEu0', 'nuEu', 'nuAs0', 'nuAs', 'nuJp0', 'nuJp',
              'mAfB', 'mAfEu', 'mAfAs', 'mEuAs', 'mChJp',
              'TAf', 'TB', 'TEuAs', 'TChJp']

# Parameters from Jouganous et al. (2017) paper (Table 4):
# N_A = 11293
# N_Af = 23721
# N_B = 2831
# N_Eu0 = 2512
# r_Eu = 0.16e-2  # it is percent in table
# N_As0 = 1019
# r_As = 0.26e-2  # it is percent in table
# N_Jp0 = 4384
# r_Jp = 1.29e-2  # it is percent in table
# m_Af_B = 16.8e-5
# m_Af_Eu = 1.14e-5
# m_Af_As = 0.56e-5
# m_Eu_As = 4.75e-5
# m_Ch_Jp = 3.3e-5
# T_Af = 357e3  # kya
# T_B = 119e3  # kya
# T_Eu_As = 46e3  # kya
# T_Ch_Jp = 9e3  # kya
# We traslated them and got:
popt = [2.100504737447977, 0.25068626582838927, 0.2224386788275923,
        2.808917103428661, 0.090232887629505, 5.547927952720224,
        0.38820508279465155, 20.730881823465822, 3.7944479999999996,
        0.2574804, 0.1264816, 1.072835,
        0.7453380000000001, 0.3633621071338058, 0.11145140260826816,
        0.05648906707542359, 0.013740583883211146]

lower_bound = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2,
               0, 0, 0, 0, 0, 1e-15, 1e-15, 1e-15, 1e-15]
upper_bound = [100, 100, 100, 100, 100, 100, 100, 100,
               10, 10, 10, 10, 10, 5, 5, 5, 5]

mu = 1.44e-8  # mutation rate
L = None  # effective length of sequence

ns = [40, 40, 40, 40]

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
        fig_title=f'Demographic model for 4_YRI_CEU_CHB_JPT_17_Jou',
        pop_labels=pop_labels,
        nref=Nanc,
        draw_scale=True,
        draw_ancestors=True,
        gen_time=1.0,
        gen_time_units='Generations',
        reverse_timeline=True)
    print('Model plot is saved to model_plot.png')
