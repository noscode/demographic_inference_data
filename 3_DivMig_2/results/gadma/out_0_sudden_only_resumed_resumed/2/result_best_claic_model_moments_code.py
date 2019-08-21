#current best params = [9963.7712579862218, 12810.297488614387, 5265.214081182894, 16447.118461460097, 4987.8652873009905, 10285.17315516716, 1200.2813260233115, 464.35249445412376, 0.00016500580383414798, 4.231983944538797e-05, 2.8233855522495605e-05, 5.8959557592557226e-05, 4.5063377121153667e-05, 0.00025035171152177173, 9.4543495015313888e-05, 0.00023623017917296158]
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

popt = [1.285687633419585, 0.52843586477988325, 1.6506920959548628, 0.50060014006273845, 1.0322570529630863, 0.12046456055093144, 0.046604090201482012, 1.6440800856435964, 0.42166519990854823, 0.28131567815717728, 0.587459545324305, 0.44900118174754478, 2.4944471876482872, 0.94200975826314814, 2.3537434695124899]
ns = [20, 20, 20]
model = generated_model(popt, ns)
N_A = 9963.771258
ll_model = moments.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9963.77125799,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)