import moments
from demographic_model import model_func

n_pop = 1
pop_labels = ["A.thaliana"]

par_labels = ['N1', 'T1', 'N2', 'T2']
popt = [0.149, 0.023, 1.256, 0.045]

lower_bound = [1e-3, 0, 1e-3, 0]
upper_bound = [100, 5, 100, 5]

mu = 7e-9  # mutation rate
L = 31650335 / (2.31+1)  # effective length of sequence

ns = [16]

# Get maximum log-likelihood
model = model_func(popt, ns)
data = moments.Spectrum.from_file("fs_data.fs")
max_ll = moments.Inference.ll_multinom(model, data)

# Get ancestral population size
theta =  moments.Inference.optimal_sfs_scaling(model, data) # mutation flux
Nanc = int(theta / (4 * mu * L))

if __name__ == "__main__":
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
        fig_title=f'Demographic model for 1_AraTha_4_Hub, Nanc: {Nanc}',
        pop_labels=pop_labels,
        nref=Nanc,
        draw_scale=False,
        draw_ancestors=True,
        gen_time=1.0,
        gen_time_units='Generations',
        reverse_timeline=True)
    print('Model plot is saved to model_plot.png')
