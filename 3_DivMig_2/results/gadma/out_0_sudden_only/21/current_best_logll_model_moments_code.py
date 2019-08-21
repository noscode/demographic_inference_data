#current best params = [10146.643392381144, 13610.593371533576, 3948.7523747264049, 16093.253100707923, 4780.1558925595518, 9256.5025559945643, 752.86099690776905, 656.1967816497978, 0.0, 0.0, 3.9287182230269284e-05, 6.8429189839676888e-05, 4.1396956547599682e-05, 0.00055436018613841723, 8.91885711705666e-05, 0.00046332742648955199]
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

popt = [10146.643392381144, 13610.593371533576, 3948.7523747264049, 16093.253100707923, 4780.1558925595518, 9256.5025559945643, 752.86099690776905, 656.1967816497978, 0.0, 0.0, 3.9287182230269284e-05, 6.8429189839676888e-05, 4.1396956547599682e-05, 0.00055436018613841723, 8.91885711705666e-05, 0.00046332742648955199]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 1.3413887573651622, 0.38916834090093505, 1.5860666900732845, 0.47110711470838218, 0.91227238388460918, 0.074198034541459609, 0.064671316047484162, 0.0, 0.0, 0.3986330279820357, 0.6943265869327524, 0.42004015561839164, 5.6248951196805521, 0.90496462634374497, 4.7012181704991729]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10146.6433924,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)