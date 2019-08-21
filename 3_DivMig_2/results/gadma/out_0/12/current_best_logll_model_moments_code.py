#current best params = [10181.573503521164, 0.83235941722307161, 31435.891849720541, 36932.481823922746, 0.9905129831065187, 5997.2570862716902, 3201.6351794020802, 19967.453777624971, 1115.1785505183107, 244.41849122735593, 0.0, 0.0, 8.3080501664264096e-05, 0.00015559504482373686, 8.3635300218185633e-05, 4.0063279730723244e-05, 0.00017223517163989316, 0.00047262717221197253]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

def generated_model((s0, nu21, nu22, t1, m1_12, m1_21, s1, nu31, nu32, nu33,
		t2, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32), ns):
	theta1 = 2.0
	sts = moments.LinearSystem_1D.steady_state_1D(sum(ns), theta=theta1)
	fs = moments.Spectrum(sts)

	fs = moments.Manips.split_1D_to_2D(fs, ns[0], sum(ns[1:]))

	before = [s0, 1 - s0]
	T = nu31
	after = nu21, nu22
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T), lambda t: before[1] * (after[1] / before[1]) ** (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, nu33],[t2, 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - t1) * before[-1])
	before[-2] *= t1
	T = nu32
	after = m1_12, m1_21, s1
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T), lambda t: after[1], lambda t: before[2] + (after[2] - before[2]) * (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, m2_12, m2_13],[m2_21, 0.0, m2_23], [m2_31, m2_32, 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.83235941722307161, 3.0875278598979667, 3.6273844913215161, 0.9905129831065187, 0.58903047590803304, 0.31445386887349353, 1.9611363381818623, 0.10952909686627914, 0.024005964416288602, 0.0, 0.0, 0.84589023440411726, 1.584202385656547, 0.85153895666051671, 0.40790722737048829, 1.753625059943156, 4.8120882936375535]
ns = [20, 20, 20]
model = generated_model(popt, ns)
N_A = 10181.573504
ll_model = moments.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10181.5735035,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)