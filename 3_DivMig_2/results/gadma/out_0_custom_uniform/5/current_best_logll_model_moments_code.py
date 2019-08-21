#current best params = [1.495654698470158, 0.49388075356037237, 0.98022956243691772, 0.49248022571803585, 0.98462990682331231, 3.3853614546916, 0.098151534612967703, 0.051689257434825164]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

import imp
file_with_model_func = imp.load_source("file_with_model_func= file_with_model_func", "/mnt/ssd1/enoskova/simulations/3_DivMig_2/demographic_model.py")
generated_model = file_with_model_func.model_func
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [1.495654698470158, 0.49388075356037237, 0.98022956243691772, 0.49248022571803585, 0.98462990682331231, 3.3853614546916, 0.098151534612967703, 0.051689257434825164]
ns = [20, 20, 20]
model = generated_model(popt, ns)
Nref = 0.0516892574348# It is also optimized parameter by GADMA
ll_model = Nref * dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
theta0 = 2.0
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=int(theta / theta0),
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)