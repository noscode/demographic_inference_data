#current best params = [9996.9882692370047, 0.48667262701538139, 14693.7077739742, 4661.0841320895515, 0.54289634907552176, 15481.018211297387, 4972.218753147964, 10142.326190192493, 1012.621416177032, 512.63940465703138, 2.3143394563933557e-05, 2.7811199599886778e-05, 4.5270475247589694e-05, 9.3721328691165044e-05, 4.4757048960255865e-05, 0.00034067378168311287, 9.3019688674096033e-05, 0.00029748422713635892]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

def generated_model(params, ns):
	Ns = params[:8]
	Ts = params[8:10]
	Ms = params[10:]
	theta1 = 2.0
	sts = moments.LinearSystem_1D.steady_state_1D(sum(ns), theta=theta1, N=Ns[0])
	fs = moments.Spectrum(sts)

	fs = moments.Manips.split_1D_to_2D(fs, ns[0], sum(ns[1:]))

	before.append((1 - Ns[1]) * before[-1])
	before[-2] *= Ns[1]
	T = Ts[0]
	after = Ns[2:4]
	growth_funcs = [lambda t: after[0], lambda t: before[1] * (after[1] / before[1]) ** (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, Ms[0]],[Ms[1], 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - Ns[4]) * before[-1])
	before[-2] *= Ns[4]
	T = Ts[1]
	after = Ns[5:8]
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: after[1], lambda t: after[2]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [9996.9882692370047, 0.48667262701538139, 14693.7077739742, 4661.0841320895515, 0.54289634907552176, 15481.018211297387, 4972.218753147964, 10142.326190192493, 1012.621416177032, 512.63940465703138, 2.3143394563933557e-05, 2.7811199599886778e-05, 4.5270475247589694e-05, 9.3721328691165044e-05, 4.4757048960255865e-05, 0.00034067378168311287, 9.3019688674096033e-05, 0.00029748422713635892]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.48667262701538139, 1.4698134456344281, 0.46624883480485441, 0.54289634907552176, 1.5485682081808563, 0.4973716702707961, 1.0145381706011125, 0.10129264823617902, 0.051279384435664377, 0.23136424396596722, 0.27802823615347699, 0.45256840999293835, 0.93693102350288249, 0.44743569342134415, 3.4057117991226877, 0.92991673648301632, 2.9739463289652166]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9996.98826924,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)