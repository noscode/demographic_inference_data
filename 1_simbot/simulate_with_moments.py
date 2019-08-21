import moments
from demographic_model import model_func as func

theta = 20000
popt = [0.01, 1.0, 0.005, 0.05]
p = [ 0.00419757,  0.99078189,  0.00207337,  0.05093703]
p = [ 29.7588,   0.0589,  29.8005,   0.3624]
p = [ 0.03946682,  1.01402808,  0.02096663,  0.0464761 ]
p = [ 0.02194382,  1.01816194,  0.01115241,  0.04853506]
p = [ 0.0202,  1.013 ,  0.0103,  0.0486]

data = func(popt, [20]) * theta
model = data
ll_model = moments.Inference.ll_multinom(model, model)
print('Maximum log composite likelihood: {0}'.format(ll_model))

model = func(p, [20]) 
ll_model = moments.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))
theta = moments.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta /2: {0}'.format(theta / 2))
print([x * theta / 2 for x in p])


model *= theta
i = 0
def func_2(p, ns):
	print 'something'
	return func(p, ns)

lower_bound = [1e-3, 1e-3, 0, 0]
upper_bound = [100, 100, 5, 5]
p0 = moments.Misc.perturb_params(popt, fold=1, upper_bound=upper_bound,
                              lower_bound=lower_bound)
# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters. If you want to fix
# theta, add a multinom=False to the call.
# The maxiter argument restricts how long the optimizer will run. For real 
# runs, you will want to set this value higher (at least 10), to encourage
# better convergence. You will also want to run optimization several times
# using multiple sets of intial parameters, to be confident you've actually
# found the true maximum likelihood parameters.
print('Beginning optimization ************************************************')
popt = moments.Inference.optimize_powell(p0, data, func_2, 
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=1, maxiter=None)
# The verbose argument controls how often progress of the optimizer should be
# printed. It's useful to keep track of optimization process.
print('Finshed optimization **************************************************')
#model.to_file('simulated_data/00.fs')
