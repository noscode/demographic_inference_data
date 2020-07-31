import moments
from demographic_model import model_func as func

theta = 20000
par_labels = ['nuB', 'nuF', 'tB', 'tF']
popt = [0.01, 1.0, 0.005, 0.05]

lower_bound = [1e-3, 1e-3, 0, 0]
upper_bound = [100, 100, 5, 5]

data = func(popt, [20]) * theta
model = data
ll_model = moments.Inference.ll_multinom(model, model)
print('Maximum log composite likelihood: {0}'.format(ll_model))

model = func(p, [20])
ll_model = moments.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))
theta = moments.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))

mu = 2.5e-8  # mutation rate
L = 20000000  # effective length of sequence
theta0 = 4 * mu * L

Nanc = theta / theta0
print('Size of the ancestral population: {0}'.format(Nanc))
