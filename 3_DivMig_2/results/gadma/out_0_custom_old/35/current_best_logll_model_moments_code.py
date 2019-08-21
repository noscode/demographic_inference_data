#current best params = [1.6053302589746288, 0.60428549607958193, 1.0456403396755352, 0.61576729345551717, 1.4059904967126731, 0.0069628483336063134, 0.1351544161066752, 0.041326740898936058]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

import imp
file_with_model_func = imp.load_source("file_with_model_func= file_with_model_func", "/mnt/ssd1/enoskova/simulations/3_DivMig_2/demographic_model.py")
generated_model = file_with_model_func.model_func
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [1.6053302589746288, 0.60428549607958193, 1.0456403396755352, 0.61576729345551717, 1.4059904967126731, 0.0069628483336063134, 0.1351544161066752, 0.041326740898936058]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll_multinom(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))

theta = moments.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.6053302589746288, 0.60428549607958193, 1.0456403396755352, 0.61576729345551717, 1.4059904967126731, 0.0069628483336063134, 0.1351544161066752, 0.041326740898936058]
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