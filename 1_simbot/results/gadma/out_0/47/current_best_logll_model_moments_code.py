#current best params = [10322.981388681661, 333.16687386803198, 27881.401013195049, 486.35774426575489, 727.53664135870872]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

def generated_model(params, ns):
	Ns = params[:3]
	Ts = params[3:5]
	Ms = params[5:]
	theta1 = 2.0
	sts = moments.LinearSystem_1D.steady_state_1D(sum(ns), theta=theta1, N=Ns[0])
	before = [Ns[0]]
	fs = moments.Spectrum(sts)

	T = Ts[0]
	after = Ns[1:2]
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	fs.integrate(Npop=list_growth_funcs, tf=T, dt_fac=0.1, theta=theta1)

	before = after
	T = Ts[1]
	after = Ns[2:3]
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	fs.integrate(Npop=list_growth_funcs, tf=T, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/1_simbot/simulated_data/00.fs')

popt = [10322.981388681661, 333.16687386803198, 27881.401013195049, 486.35774426575489, 727.53664135870872]
ns = [20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.032274287952637724, 2.7009058685085705, 0.047114077411687286, 0.07047737605692049]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0'],
	nref=10322.9813887,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)