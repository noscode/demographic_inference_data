#current best params = [5478.4248493175191, 3596.7223430347676, 9725.3644041763546, 7759.2240293939094, 5052.1305913027772, 5197.8683558812763, 5475.5417105266342, 27392.124246587595, 0.00071772826014466747, 0.00018254949650226883, 2.5938824407668029e-17, 0.0004618919277153057, 0.0001247936450581875, 0.0, 3.7890755505532726e-07, 0.0018253421877724141]
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

popt = [0.65652490304450073, 1.7752117938403231, 1.4163238965230167, 0.92218671064405533, 0.9487888396477332, 0.9994737285131795, 5.0, 3.9320203354339749, 1.0000836978684311, 1.421039001970523e-13, 2.530440214494702, 0.68367260612368475, 0.0, 0.0020758165652092508, 9.9999999999999982]
ns = [20, 20, 20]
model = generated_model(popt, ns)
N_A = 5478.424849
ll_model = moments.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=5478.42484932,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)