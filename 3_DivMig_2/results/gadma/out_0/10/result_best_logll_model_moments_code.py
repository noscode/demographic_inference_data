#current best params = [9946.2946761870062, 0.38980219562901031, 12638.92408523819, 5053.8278596744503, 0.76314877856248642, 20775.198406385371, 5875.0628557831023, 10570.420137094325, 1258.5161692946201, 449.48288438495638, 0.00020677403051294664, 3.0423300031225863e-05, 1.6187493624783026e-05, 5.2938848903833364e-05, 4.8691633470202338e-05, 0.00024180729545905956, 9.5511813893264647e-05, 0.00018566320755764254]
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
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T), lambda t: before[1] + (after[1] - before[1]) * (t / T), lambda t: after[2]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, m2_12, m2_13],[m2_21, 0.0, m2_23], [m2_31, m2_32, 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.38980219562901031, 1.2707168344305908, 0.50811161585370168, 0.76314877856248642, 2.0887374728726331, 0.5906785438248604, 1.0627495445516584, 0.12653115660324299, 0.04519098810344812, 2.0566354388646508, 0.30259910713262178, 0.16100558166099052, 0.52654539101566611, 0.48430133475952253, 2.4050866154876225, 0.94998864603954225, 1.846660972894383]
ns = [20, 20, 20]
model = generated_model(popt, ns)
N_A = 9946.294676
ll_model = moments.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9946.29467619,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)