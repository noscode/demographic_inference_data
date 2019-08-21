#current best params = [9912.4421864740798, 9495.8443096234496, 5555.6526411173381, 18771.380061731641, 4822.3799303490187, 8841.4665586861647, 1811.6241829731266, 500.04579558926002, 0.0004417740831135985, 0.00012219149657863513, 0.0, 0.0, 4.6476516436011135e-06, 0.00042862942627585449, 8.0214671578383414e-05, 0.00037827047274219787]
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

popt = [0.95797222631783974, 0.56047264000169872, 1.8937189956422578, 0.48649766017595009, 0.89195643135762204, 0.18276264808335121, 0.050446276122708925, 4.3790600583461403, 1.2112161455144661, 0.0, 0.0, 0.04606957822006727, 4.2487644073809614, 0.79512329452773112, 3.7495841919072554]
ns = [20, 20, 20]
model = generated_model(popt, ns)
N_A = 9912.442186
ll_model = moments.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9912.44218647,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)