#current best params = [10022.49, 14707.07007694302, 4446.3390843136176, 15166.561393827767, 5058.0470727455577, 9917.3650173477236, 968.03672998745299, 543.4223675150871, 0.0, 8.9648302589830816e-05, 4.7279542926963532e-05, 9.0618210262992068e-05, 3.4995573579229827e-05, 0.00035187383071244279, 8.129780817389596e-05, 0.00037068468589426557]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

def generated_model((nu21, nu22, t1, m1_12, m1_21, nu31, nu32, nu33, t2, m2_12,
		m2_13, m2_21, m2_23, m2_31, m2_32), ns):
	theta1 = 2.0
	sts = moments.LinearSystem_1D.steady_state_1D(sum(ns), theta=theta1)
	fs = moments.Spectrum(sts)

	fs = moments.Manips.split_1D_to_2D(fs, ns[0], sum(ns[1:]))

	T = nu31
	after = nu21, nu22
	growth_funcs = [lambda t: after[0], lambda t: after[1]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, nu33],[t2, 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	T = nu32
	after = t1, m1_12, m1_21
	growth_funcs = [lambda t: after[0], lambda t: after[1], lambda t: after[2]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, m2_12, m2_13],[m2_21, 0.0, m2_23], [m2_31, m2_32, 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [1.4674068097791089, 0.44363617068349459, 1.5132528337596514, 0.50466970510776843, 0.98951109129046011, 0.096586450072532179, 0.05422029530736245, 0.0, 0.89849921622355344, 0.47385874619006274, 0.90822010617873539, 0.35074278624209515, 3.5266519495771509, 0.81480646944479052, 3.7151835575284178]
ns = [20, 20, 20]
model = generated_model(popt, ns)
N_A = 10022.490000
ll_model = moments.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10022.49,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)