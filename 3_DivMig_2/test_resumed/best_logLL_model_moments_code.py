#current best params = [9989.4660170128627, 13260.272953736603, 5179.8271580520213, 16112.333171372218, 4961.5111076850453, 10095.153270706542, 1114.7427994024242, 482.66679968862951, 0.00011989613988540882, 1.9311681448141275e-05, 3.5280234702241428e-05, 7.0423934210191519e-05, 4.8628650245840641e-05, 0.00028585742083221362, 9.8696522105543903e-05, 0.00027235644764340468]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

def generated_model(params, ns):
	Ns = params[:6]
	Ts = params[6:8]
	Ms = params[8:]
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
	T = Ts[1]
	after = Ns[3:6]
	growth_funcs = [lambda t: after[0], lambda t: after[1], lambda t: after[2]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [9989.4660170128627, 13260.272953736603, 5179.8271580520213, 16112.333171372218, 4961.5111076850453, 10095.153270706542, 1114.7427994024242, 482.66679968862951, 0.00011989613988540882, 1.9311681448141275e-05, 3.5280234702241428e-05, 7.0423934210191519e-05, 4.8628650245840641e-05, 0.00028585742083221362, 9.8696522105543903e-05, 0.00027235644764340468]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 1.3274256032457885, 0.51852893330137562, 1.6129323773594726, 0.49667430663813195, 1.0105798701866231, 0.11159183058473073, 0.04831757762292891, 1.1976984149563119, 0.19291338555758503, 0.35243070563027867, 0.70349749757705782, 0.48577424908402927, 2.855562991114343, 0.98592555357068967, 2.7206954782481341]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9989.46601701,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)