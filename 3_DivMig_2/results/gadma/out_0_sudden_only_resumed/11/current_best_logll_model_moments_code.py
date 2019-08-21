#current best params = [9981.4360801878574, 12747.162665906686, 5325.3181404140105, 16557.061570249938, 4975.2668538190892, 10173.555653822637, 1179.2279579301801, 465.17695506418397, 0.00017008931620344906, 2.1763793223883738e-05, 2.7129241976365652e-05, 5.5405711690783118e-05, 4.9221887997828585e-05, 0.00025142413536569298, 0.00010075446750901771, 0.00024194877729613324]
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

popt = [1.2770870407323969, 0.53352224045037255, 1.6587855131501601, 0.49845200769200848, 1.019247688618286, 0.11814211386584225, 0.046604211190363003, 1.6977356376075876, 0.21723391092662114, 0.27078879469104306, 0.55302856971886882, 0.49130512879649191, 2.509573936169164, 1.0056742772346245, 2.4149962552609612]
ns = [20, 20, 20]
model = generated_model(popt, ns)
N_A = 9981.436080
ll_model = moments.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9981.43608019,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)