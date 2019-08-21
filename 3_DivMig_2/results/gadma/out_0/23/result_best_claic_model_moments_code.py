#current best params = [9982.0799086782972, 0.58467453223565291, 15254.678507303075, 6844.9324745553113, 0.048366945598932842, 14310.694221270651, 4770.8066039537016, 12726.376303673798, 1076.8212536479452, 428.84890656639442, 0.0, 9.5460547671437871e-06, 5.3643712008808427e-05, 0.0001153814109595566, 5.2301200522052453e-05, 0.00022528115111516839, 0.00010607170081641228, 0.00023346863339923779]
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
	growth_funcs = [lambda t: after[0], lambda t: before[1] * (after[1] / before[1]) ** (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, nu33],[t2, 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - t1) * before[-1])
	before[-2] *= t1
	T = nu32
	after = m1_12, m1_21, s1
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: after[1], lambda t: before[2] + (after[2] - before[2]) * (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, m2_12, m2_13],[m2_21, 0.0, m2_23], [m2_31, m2_32, 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.58467453223565291, 1.5282064105738971, 0.68572206766291377, 0.048366945598932842, 1.4336385154389628, 0.47793712809351702, 1.2749223027767633, 0.10787543913686466, 0.042961878735668951, 0.0, 0.095289481498248674, 0.53547581987005133, 1.1517464641743438, 0.52207476293093469, 2.2487744523506419, 1.0588161935988443, 2.3305025547611105]
ns = [20, 20, 20]
model = generated_model(popt, ns)
N_A = 9982.079909
ll_model = moments.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9982.07990868,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)