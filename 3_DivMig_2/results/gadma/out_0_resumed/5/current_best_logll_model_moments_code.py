#current best params = [10060.2253535073, 0.78119533145150966, 34778.113098094385, 10546.024389682447, 0.50945947763997879, 3348.5808421657448, 4285.5350216377146, 12163.809500597694, 1112.3323248553168, 341.99366805514808, 5.1461772907270248e-06, 0.0, 6.7054119021242473e-05, 0.00014151304977544829, 6.7947610583648195e-05, 0.00011912166643842915, 0.00013839309354127111, 0.00017893863068120855]
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
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T), lambda t: before[1] + (after[1] - before[1]) * (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, Ms[0]],[Ms[1], 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - Ns[4]) * before[-1])
	before[-2] *= Ns[4]
	T = Ts[1]
	after = Ns[5:8]
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: after[1], lambda t: before[2] + (after[2] - before[2]) * (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [10060.2253535073, 0.78119533145150966, 34778.113098094385, 10546.024389682447, 0.50945947763997879, 3348.5808421657448, 4285.5350216377146, 12163.809500597694, 1112.3323248553168, 341.99366805514808, 5.1461772907270248e-06, 0.0, 6.7054119021242473e-05, 0.00014151304977544829, 6.7947610583648195e-05, 0.00011912166643842915, 0.00013839309354127111, 0.00017893863068120855]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.78119533145150966, 3.4569914565551634, 1.048289080920616, 0.50945947763997879, 0.33285346247222264, 0.42598797452818959, 1.2090991079396667, 0.11056733679107138, 0.033994632926976999, 0.051771703253815519, 0.0, 0.67457954823459965, 1.4236531712031053, 0.6835682747038585, 1.1983908088559245, 1.3922657083942029, 1.8001629491009734]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10060.2253535,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)