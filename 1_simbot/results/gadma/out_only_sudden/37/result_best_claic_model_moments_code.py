#current best params = [10121.467746376455, 517.01279597233838, 10685.044394328719, 284.82790550735621, 442.93421227022503]
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
	fs = moments.Spectrum(sts)

	T = Ts[0]
	after = Ns[1:2]
	growth_funcs = [lambda t: after[0]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	fs.integrate(Npop=list_growth_funcs, tf=T, dt_fac=0.1, theta=theta1)

	T = Ts[1]
	after = Ns[2:3]
	growth_funcs = [lambda t: after[0]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	fs.integrate(Npop=list_growth_funcs, tf=T, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/1_simbot/SimBot.fs')

popt = [10121.467746376455, 517.01279597233838, 10685.044394328719, 284.82790550735621, 442.93421227022503]
ns = [20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.051080812479734672, 1.0556813163934675, 0.028140968547701619, 0.043761855826572001]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop_1'],
	nref=10121.4677464,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)