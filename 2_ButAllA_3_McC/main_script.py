import moments
from demographic_model import model_func

n_pop = 2
pop_labels = ["WY", "CO"]

par_labels = ['nuW', 'nuC', 'T']
popt = [1.320, 0.173, 0.117]

lower_bound = [1e-2, 1e-2, 1e-15]
upper_bound = [100, 100, 5]

mu = None  # mutation rate
L = None  # effective length of sequence
Nanc = None

ns = [12, 12]

# Get maximum log-likelihood
model = model_func(popt, ns)
data = moments.Spectrum.from_file("fs_data.fs")
max_ll = moments.Inference.ll_multinom(model, data)

if __name__ == "__main__":
    print('Maximum log composite likelihood: {0}'.format(max_ll))

    theta = round(moments.Inference.optimal_sfs_scaling(model, data), 5)
    print('Optimal value of theta: {0}'.format(theta))

    # Draw model
    model = moments.ModelPlot.generate_model(model_func, popt,
                                             [4 for _ in range(n_pop)])
    moments.ModelPlot.plot_model(model,
        save_file='model_plot.png',
        fig_title=f'Demographic model for 2_ButAllA_3_McC',
        pop_labels=pop_labels,
        nref=None,
        draw_scale=False,
        draw_ancestors=True,
        gen_time=1.0,
        gen_time_units='Generations',
        reverse_timeline=True)
    print('Model plot is saved to model_plot.png')
