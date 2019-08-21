#current best params = [9997.7715473011449, 14933.714582964159, 4803.6442891832949, 0.99791248490025131, 15028.771811051089, 5163.6246455597984, 10079.640575784822, 1005.5757848057365, 508.09677841654081, 0.0, 4.4016925212442216e-05, 4.7887023825785343e-05, 0.00010086760240868096, 4.2704766759326242e-05, 0.00031811426407407818, 8.9486321161522117e-05, 0.00031012211652042154]
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
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: before[1] * (after[1] / before[1]) ** (t / T), lambda t: after[2]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [9997.7715473011449, 14933.714582964159, 4803.6442891832949, 0.99791248490025131, 15028.771811051089, 5163.6246455597984, 10079.640575784822, 1005.5757848057365, 508.09677841654081, 0.0, 4.4016925212442216e-05, 4.7887023825785343e-05, 0.00010086760240868096, 4.2704766759326242e-05, 0.00031811426407407818, 8.9486321161522117e-05, 0.00031012211652042154]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 1.4937043232394573, 0.48047149971935676, 0.99791248490025131, 1.5032121648256747, 0.51647755913703552, 1.0081887276675949, 0.10057999225609304, 0.050821003061797237, 0.0, 0.4400711624886372, 0.47876352429036872, 1.0084512454059948, 0.4269525020405236, 3.1804337381504615, 0.89466379558131814, 3.1005300727366807]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9997.7715473,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)