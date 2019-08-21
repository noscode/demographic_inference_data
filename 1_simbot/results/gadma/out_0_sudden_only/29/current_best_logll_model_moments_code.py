#current best params = [67249.884359338394, 2783.2480207563062, 6724988.4359338395, 8295.594770911106, 147.4286048376965]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

def generated_model((nu21, t1, nu31, t2), ns):
	theta1 = 2.0
	sts = moments.LinearSystem_1D.steady_state_1D(sum(ns), theta=theta1)
	fs = moments.Spectrum(sts)

	T = nu31
	after = nu21
	growth_funcs = [lambda t: after[0]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	fs.integrate(Npop=list_growth_funcs, tf=T, dt_fac=0.1, theta=theta1)

	T = t2
	after = t1
	growth_funcs = [lambda t: after[0]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	fs.integrate(Npop=list_growth_funcs, tf=T, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/1_simbot/simulated_data/00.fs')

popt = [0.04138665883623667, 100.0, 0.1233547812005891, 0.0021922506817994908]
ns = [20]
model = generated_model(popt, ns)
N_A = 67249.884359
ll_model = moments.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0'],
	nref=67249.8843593,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)