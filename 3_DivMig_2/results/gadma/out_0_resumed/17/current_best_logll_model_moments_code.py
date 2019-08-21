#current best params = [10126.092521942077, 0.79617241371674641, 7365.7859840467172, 3648.4156753471652, 0.57812963138303286, 15956.680256129115, 4863.4845393689393, 12431.054161907363, 5564.9262054603114, 1217.3504614242142, 0.00067645276931379032, 0.00098754776122489006, 2.9595426134274518e-05, 6.8460559673336394e-05, 3.552270116437394e-05, 0.00076506501994233532, 8.7532955250038197e-05, 0.00079792922809833759]
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
	growth_funcs = [lambda t: after[0], lambda t: before[1] + (after[1] - before[1]) * (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, Ms[0]],[Ms[1], 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - Ns[4]) * before[-1])
	before[-2] *= Ns[4]
	T = Ts[1]
	after = Ns[5:8]
	growth_funcs = [lambda t: after[0], lambda t: before[1] + (after[1] - before[1]) * (t / T), lambda t: before[2] * (after[2] / before[2]) ** (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [10126.092521942077, 0.79617241371674641, 7365.7859840467172, 3648.4156753471652, 0.57812963138303286, 15956.680256129115, 4863.4845393689393, 12431.054161907363, 5564.9262054603114, 1217.3504614242142, 0.00067645276931379032, 0.00098754776122489006, 2.9595426134274518e-05, 6.8460559673336394e-05, 3.552270116437394e-05, 0.00076506501994233532, 8.7532955250038197e-05, 0.00079792922809833759]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.79617241371674641, 0.72740654582070097, 0.36029847322068886, 0.57812963138303286, 1.5757983863521714, 0.48029232686056622, 1.2276259707256969, 0.54956304155840541, 0.12021917228055698, 6.8498233287953809, 10.0, 0.29968602326196631, 0.69323796135614102, 0.35970615861975008, 7.7471191772375478, 0.88636680358090225, 8.07990518968559]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10126.0925219,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)