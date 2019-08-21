#current best params = [9989.5014354704326, 15083.152331413134, 4765.1085666310355, 14863.191908376666, 5075.0491814587249, 9989.7614998671161, 1012.5322772394468, 512.47583389916963, 4.6552426410199958e-06, 5.1999221783132496e-05, 5.2200762229293804e-05, 9.4397442791185497e-05, 4.1125518414654727e-05, 0.00030043872576316337, 8.9808258547042707e-05, 0.0003280890990819378]
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

popt = [9989.5014354704326, 15083.152331413134, 4765.1085666310355, 14863.191908376666, 5075.0491814587249, 9989.7614998671161, 1012.5322772394468, 512.47583389916963, 4.6552426410199958e-06, 5.1999221783132496e-05, 5.2200762229293804e-05, 9.4397442791185497e-05, 4.1125518414654727e-05, 0.00030043872576316337, 8.9808258547042707e-05, 0.0003280890990819378]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 1.5099004118319972, 0.47701165042243504, 1.4878812525719125, 0.50803828541816776, 1.0000260337713911, 0.10135964079690468, 0.051301442540414001, 0.046503553044932415, 0.51944630064594743, 0.52145958922218116, 0.94298339026728561, 0.41082342523765913, 3.0012330822820283, 0.89713972767278294, 3.2774465262412189]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9989.50143547,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)