#current best params = [9998.4978265828595, 9997.6928325898934, 1005.1705787227135, 503.10153897513879, 0.00050082622903874078, 0.00024989272447790999]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

def generated_model(params, ns):
	Ns = params[:3]
	Ts = params[3:4]
	Ms = params[4:]
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

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/2_DivMigr/simulated_data/00.fs')

popt = [9998.4978265828595, 9997.6928325898934, 1005.1705787227135, 503.10153897513879, 0.00050082622903874078, 0.00024989272447790999]
ns = [20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.99991948850648094, 0.10053215954603513, 0.050317712490525343, 5.0075099625395385, 2.4985518625712522]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1'],
	nref=9998.49782658,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)