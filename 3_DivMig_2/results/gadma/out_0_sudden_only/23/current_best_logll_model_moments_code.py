#current best params = [10191.272896673674, 12526.969362163436, 2712.8963465154347, 15875.645189226432, 4439.7486717685251, 8345.9850404495428, 491.02748168904128, 874.65671419047089, 0.0, 7.0549241298098844e-05, 3.6550655965642258e-05, 4.7063913869185431e-05, 3.5904173804096048e-05, 0.00072376613894435931, 8.2721456526511618e-05, 0.00060487436724135021]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

def generated_model(params, ns):
	Ns = params[:6]
	Ts = params[6:8]
	Ms = params[8:]
	theta1 = 2.0
	sts = moments.LinearSystem_1D.steady_state_1D(sum(ns), theta=theta1, N=Ns[0])
	fs = moments.Spectrum(sts)

	fs = moments.Manips.split_1D_to_2D(fs, ns[0], sum(ns[1:]))

	T = Ts[0]
	after = Ns[1:3]
	growth_funcs = [lambda t: after[0], lambda t: after[1]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, Ms[0]],[Ms[1], 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	T = Ts[1]
	after = Ns[3:6]
	growth_funcs = [lambda t: after[0], lambda t: after[1], lambda t: after[2]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [10191.272896673674, 12526.969362163436, 2712.8963465154347, 15875.645189226432, 4439.7486717685251, 8345.9850404495428, 491.02748168904128, 874.65671419047089, 0.0, 7.0549241298098844e-05, 3.6550655965642258e-05, 4.7063913869185431e-05, 3.5904173804096048e-05, 0.00072376613894435931, 8.2721456526511618e-05, 0.00060487436724135021]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 1.2291859406740162, 0.26619798861444444, 1.5577686271562876, 0.43564221238915241, 0.81893450652014088, 0.048181172917988249, 0.085824089204396617, 0.0, 0.71898657072220584, 0.37249770949829392, 0.47964118982641374, 0.36590923336714498, 7.376098235353802, 0.84303693787200751, 6.1644397447594113]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10191.2728967,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)