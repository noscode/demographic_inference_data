#current best params = [10547.124587666081, 0.49922884433434311, 2324.0720055145621, 983.10038939237995, 0.38259572484471593, 245948.91429747228, 12483.945334267592, 13623.867238614041, 0.00055736129051751568, 1068.6454785185476, 2.3685334416298605e-08, 0.0, 2.3564818400192407e-12, 4.9316767303062735e-05, 2.6074236031866341e-07, 0.00094812571112453769, 1.6269337277445202e-06, 0.00047524551485529274]
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
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T), lambda t: after[1]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, nu33],[t2, 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - t1) * before[-1])
	before[-2] *= t1
	T = nu32
	after = m1_12, m1_21, s1
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T), lambda t: before[1] * (after[1] / before[1]) ** (t / T), lambda t: before[2] + (after[2] - before[2]) * (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, m2_12, m2_13],[m2_21, 0.0, m2_23], [m2_31, m2_32, 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.49922884433434311, 0.22035124229331246, 0.093210275579946017, 0.38259572484471593, 23.319048926859889, 1.1836349547692315, 1.2917138813877229, 5.2844856992520962e-08, 0.10132102542604199, 0.00024981217298923668, 0.0, 2.4854107555255544e-08, 0.52015008900633963, 0.0027500821595630638, 10.0, 0.017159472722397463, 5.012473654908284]
ns = [20, 20, 20]
model = generated_model(popt, ns)
N_A = 10547.124588
ll_model = moments.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10547.1245877,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)