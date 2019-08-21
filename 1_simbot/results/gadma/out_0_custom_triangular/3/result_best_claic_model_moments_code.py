#current best params = [4.0802262409526007, 0.18114923735150498, 4.9673007106637685, 0.57503449752124958]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

import imp
file_with_model_func = imp.load_source("file_with_model_func= file_with_model_func", "/mnt/ssd1/enoskova/simulations/1_simbot/demographic_model.py")
generated_model = file_with_model_func.model_func
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/1_simbot/simulated_data/00.fs')

popt = [4.0802262409526007, 0.18114923735150498, 4.9673007106637685, 0.57503449752124958]
ns = [20]
model = generated_model(popt, ns)
Nref = 0.575034497521# It is also optimized parameter by GADMA
ll_model = Nref * dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
theta0 = 2.0
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0'],
	nref=int(theta / theta0),
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)