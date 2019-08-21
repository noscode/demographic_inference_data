#current best params = [9967.3146445165075, 15330.245064992352, 5041.7803921540635, 0.80046644362610975, 14072.792665630384, 5768.805130313388, 10599.895470591055, 1070.4862053025111, 462.62585227749247, 2.0353352299054581e-10, 4.5591126143998023e-05, 4.7862183191897412e-05, 0.00011366669235346279, 4.388469326274478e-05, 0.00026237246902314663, 8.9994004175019657e-05, 0.00020443544383591698]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

def generated_model((nu21, nu22, t1, m1_12, m1_21, s1, nu31, nu32, nu33, t2,
		m2_12, m2_13, m2_21, m2_23, m2_31, m2_32), ns):
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
	before.append((1 - t1) * before[-1])
	before[-2] *= t1
	T = nu32
	after = m1_12, m1_21, s1
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: before[1] + (after[1] - before[1]) * (t / T), lambda t: after[2]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, m2_12, m2_13],[m2_21, 0.0, m2_23], [m2_31, m2_32, 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [1.5380516831006481, 0.50583136701998122, 0.80046644362610975, 1.4118940925951902, 0.57877225070717331, 1.0634655219219511, 0.10739966013729046, 0.046414291991073532, 2.0286826643537044e-06, 0.45442109927507091, 0.47705743944713092, 1.1329516872884222, 0.43741254582787098, 2.6151489528123633, 0.8969985557323531, 2.0376723932039673]
ns = [20, 20, 20]
model = generated_model(popt, ns)
N_A = 9967.314645
ll_model = moments.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9967.31464452,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)