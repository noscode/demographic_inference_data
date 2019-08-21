#current best params = [3121.3506828137083, 118.38763233475149, 2958.1586039585686, 0.021729657495044304, 13611.963344370251, 1935.8475687209079, 7135.4281660068946, 12455.067551026405, 6808.3795050577228, 0.0, 0.0, 0.00020079560134750286, 3.7490515077963783e-05, 0.0, 0.0029267528410053629, 0.00019315565445954418, 0.0]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

def generated_model(params, ns):
	Ns = params[:7]
	Ts = params[7:9]
	Ms = params[9:]
	theta1 = 2.0
	sts = moments.LinearSystem_1D.steady_state_1D(sum(ns), theta=theta1, N=Ns[0])
	fs = moments.Spectrum(sts)

	fs = moments.Manips.split_1D_to_2D(fs, ns[0], sum(ns[1:]))

	T = Ts[0]
	after = Ns[1:3]
	growth_funcs = [lambda t: after[0], lambda t: after[1]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, Ms[0]],[Ms[1], 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - Ns[3]) * before[-1])
	before[-2] *= Ns[3]
	T = Ts[1]
	after = Ns[4:7]
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: after[1], lambda t: before[2] * (after[2] / before[2]) ** (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [3121.3506828137083, 118.38763233475149, 2958.1586039585686, 0.021729657495044304, 13611.963344370251, 1935.8475687209079, 7135.4281660068946, 12455.067551026405, 6808.3795050577228, 0.0, 0.0, 0.00020079560134750286, 3.7490515077963783e-05, 0.0, 0.0029267528410053629, 0.00019315565445954418, 0.0]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.037928334354290563, 0.94771748020698787, 0.021729657495044304, 4.3609208729151483, 0.62019547479229686, 2.2860065693018314, 3.9902813931175833, 2.1812286400707728, 0.0, 0.0, 0.62675348737201719, 0.11702104483763988, 0.0, 9.1354219786990498, 0.60290653393662696, 0.0]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=3121.35068281,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)