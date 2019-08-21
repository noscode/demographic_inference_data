#current best params = [10071.618384336105, 14450.439230458076, 4495.1567994643083, 15481.462603724602, 4980.323396077536, 9565.2725303732823, 873.16289249195586, 571.83526452306808, 0.0, 0.0, 4.5917928887840297e-05, 7.8780344253524492e-05, 4.5055889521805194e-05, 0.00039536181315565609, 9.4259106122512016e-05, 0.00041355825188028639]
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

popt = [10071.618384336105, 14450.439230458076, 4495.1567994643083, 15481.462603724602, 4980.323396077536, 9565.2725303732823, 873.16289249195586, 571.83526452306808, 0.0, 0.0, 4.5917928887840297e-05, 7.8780344253524492e-05, 4.5055889521805194e-05, 0.00039536181315565609, 9.4259106122512016e-05, 0.00041355825188028639]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 1.4347683439765884, 0.44631921384704226, 1.5371375297342642, 0.49449087584803542, 0.94972547264595353, 0.086695390866868366, 0.056776899471530362, 0.0, 0.0, 0.46246785675741026, 0.79344556350812456, 0.4537857252304297, 3.9819333058429622, 0.94934174611458, 4.1652008926313941]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10071.6183843,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)