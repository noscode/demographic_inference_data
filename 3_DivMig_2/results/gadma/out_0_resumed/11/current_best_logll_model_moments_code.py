#current best params = [10052.490500909966, 0.66821027399324717, 38765.110144356062, 8442.6961524362887, 0.34552724848305344, 3288.2640203272285, 6736.218647503435, 13249.070838079351, 1095.9241892187292, 393.57346432158215, 1.7854070313889235e-07, 3.5413631090418735e-05, 6.2531008062886799e-05, 0.00013808924249704535, 4.8537987572488969e-05, 0.0002205692867730144, 0.00010525380506853657, 0.00022065820264739728]
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
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T), lambda t: before[1] * (after[1] / before[1]) ** (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, Ms[0]],[Ms[1], 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - Ns[4]) * before[-1])
	before[-2] *= Ns[4]
	T = Ts[1]
	after = Ns[5:8]
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: before[1] * (after[1] / before[1]) ** (t / T), lambda t: before[2] + (after[2] - before[2]) * (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [10052.490500909966, 0.66821027399324717, 38765.110144356062, 8442.6961524362887, 0.34552724848305344, 3288.2640203272285, 6736.218647503435, 13249.070838079351, 1095.9241892187292, 393.57346432158215, 1.7854070313889235e-07, 3.5413631090418735e-05, 6.2531008062886799e-05, 0.00013808924249704535, 4.8537987572488969e-05, 0.0002205692867730144, 0.00010525380506853657, 0.00022065820264739728]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.66821027399324717, 3.8562692639049985, 0.83986114203958151, 0.34552724848305344, 0.3271093884674221, 0.67010445291081477, 1.3179888940835136, 0.10902016660642699, 0.039151836481313297, 0.0017947787223295015, 0.35599519013916414, 0.62859236456449397, 1.3881407984794012, 0.48792765900573132, 2.2172706600782135, 1.058062875636093, 2.2181644860608274]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10052.4905009,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)