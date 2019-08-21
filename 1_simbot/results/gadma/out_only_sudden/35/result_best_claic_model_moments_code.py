#current best params = [9991.657474842983, 99.917827582272039, 9952.3123753309519, 50.65773536258834, 495.25124090588309]
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

popt = [9991.657474842983, 99.917827582272039, 9952.3123753309519, 50.65773536258834, 495.25124090588309]
ns = [20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.010000125387989466, 0.99606220493335618, 0.0050700031991823677, 0.049566475047090809]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop_1'],
	nref=9991.65747484,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)