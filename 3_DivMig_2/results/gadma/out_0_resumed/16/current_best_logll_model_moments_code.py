#current best params = [10257.726099746587, 0.16073514000016592, 18250.90093401623, 4943.0325055773274, 0.094575567495723808, 12643.445908071872, 15908.594382433794, 12741.4892284874, 5523.6982647987452, 354.37713391815419, 0.00031308341404593764, 0.00028856480004563624, 0.0, 0.0, 0.0, 1.1966483754996417e-05, 5.0788489746981004e-06, 0.0]
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
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T), lambda t: after[1]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, Ms[0]],[Ms[1], 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - Ns[4]) * before[-1])
	before[-2] *= Ns[4]
	T = Ts[1]
	after = Ns[5:8]
	growth_funcs = [lambda t: after[0], lambda t: before[1] + (after[1] - before[1]) * (t / T), lambda t: after[2]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [10257.726099746587, 0.16073514000016592, 18250.90093401623, 4943.0325055773274, 0.094575567495723808, 12643.445908071872, 15908.594382433794, 12741.4892284874, 5523.6982647987452, 354.37713391815419, 0.00031308341404593764, 0.00028856480004563624, 0.0, 0.0, 0.0, 1.1966483754996417e-05, 5.0788489746981004e-06, 0.0]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.16073514000016592, 1.7792345746555964, 0.48188384613812635, 0.094575567495723808, 1.2325778428012639, 1.5508889814114659, 1.2421358403011145, 0.53849149519942885, 0.034547338315741236, 3.2115239076567819, 2.9600186808962783, 0.0, 0.0, 0.0, 0.1227489127358203, 0.052097441684431901, 0.0]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10257.7260997,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)