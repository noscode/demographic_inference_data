#current best params = [9812.4441703836619, 0.13593999490561492, 32327.677077592605, 7169.326347791417, 0.07867217663034938, 9327.8912333850367, 226954.69498221355, 12090.464918790434, 1724.052437225318, 275.50892419265335, 0.00030622815919306327, 0.0, 0.0, 0.0, 6.7137815969849107e-05, 2.2945988573531492e-06, 0.00014849287969761099, 0.0]
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
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: after[1]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, nu33],[t2, 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - t1) * before[-1])
	before[-2] *= t1
	T = nu32
	after = m1_12, m1_21, s1
	growth_funcs = [lambda t: after[0], lambda t: before[1] * (after[1] / before[1]) ** (t / T), lambda t: before[2] * (after[2] / before[2]) ** (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, m2_12, m2_13],[m2_21, 0.0, m2_23], [m2_31, m2_32, 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.13593999490561492, 3.2945590839808681, 0.73063614154668866, 0.07867217663034938, 0.95061852800537472, 23.129272487197216, 1.2321563016156964, 0.17570061111062693, 0.028077502343830517, 3.0048467154812935, 0.0, 0.0, 0.0, 0.65878607092563701, 0.022515623181203922, 1.4570780917323054, 0.0]
ns = [20, 20, 20]
model = generated_model(popt, ns)
N_A = 9812.444170
ll_model = moments.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9812.44417038,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)