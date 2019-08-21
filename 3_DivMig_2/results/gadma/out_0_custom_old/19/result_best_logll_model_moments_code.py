#current best params = [1.4917100545201296, 0.49794893429435011, 1.0186914173723451, 0.51794000356066772, 0.98429484387328781, 2.8297116204349799, 0.099684906002889803, 0.049264366096534798]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

import imp
file_with_model_func = imp.load_source("file_with_model_func= file_with_model_func", "/mnt/ssd1/enoskova/simulations/3_DivMig_2/demographic_model.py")
generated_model = file_with_model_func.model_func
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [1.4917100545201296, 0.49794893429435011, 1.0186914173723451, 0.51794000356066772, 0.98429484387328781, 2.8297116204349799, 0.099684906002889803, 0.049264366096534798]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll_multinom(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))

theta = moments.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.4917100545201296, 0.49794893429435011, 1.0186914173723451, 0.51794000356066772, 0.98429484387328781, 2.8297116204349799, 0.099684906002889803, 0.049264366096534798]
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