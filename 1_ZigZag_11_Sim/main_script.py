import moments
from demographic_model import model_func

n_pop = 1
pop_labels = ["Pop 1"]

par_labels = ['nu1', 'nu2', 'nu3', 'nu4', 'nu5',
              't1', 't2', 't3', 't4', 't5', 't6']
popt = [0.1, 1.0, 0.1, 1.0, 0.1,
        1.78870737842, 0.44717719396, 0.11179429849, 0.02794857462, 0.00698693404, 0.00232902459]

lower_bound = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 0, 0, 0, 0, 0, 0]
upper_bound = [100, 100, 100, 100, 100, 5, 5, 5, 5, 5, 5]

mu = 1.25e-8  # mutation rate
L = 10000000  # effective length of sequence
Nanc = 7156
theta = 4 * mu * L * Nanc  # mutation flux

ns_per_pop = 20
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

    print('Model plot is saved to model_plot.png')
