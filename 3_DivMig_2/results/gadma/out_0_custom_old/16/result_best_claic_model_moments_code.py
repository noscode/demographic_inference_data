#current best params = [2.7823844957590107, 1.2021700571689327, 0.60893156322233677, 2.2943930972662749e-09, 1.6456129971980544, 6.8248980573684248, 0.9643881741549486, 5.1345344925699674]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

import imp
file_with_model_func = imp.load_source("file_with_model_func= file_with_model_func", "/mnt/ssd1/enoskova/simulations/3_DivMig_2/demographic_model.py")
generated_model = file_with_model_func.model_func
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [2.7823844957590107, 1.2021700571689327, 0.60893156322233677, 2.2943930972662749e-09, 1.6456129971980544, 6.8248980573684248, 0.9643881741549486, 5.1345344925699674]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll_multinom(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))

theta = moments.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [2.7823844957590107, 1.2021700571689327, 0.60893156322233677, 2.2943930972662749e-09, 1.6456129971980544, 6.8248980573684248, 0.9643881741549486, 5.1345344925699674]
theta0 = 2.0
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=int(theta / theta0),
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)