#current best params = [9989.8597311885387, 13254.579666842581, 5179.6582631480323, 16114.516866452037, 4960.964704011526, 10093.046957215538, 1114.7911879070746, 482.75216155154862, 0.00012013273452295402, 1.9328936602124865e-05, 3.5271213084713667e-05, 7.0370582534391241e-05, 4.8628468893893049e-05, 0.00028604144568468854, 9.8687003706937357e-05, 0.00027269742182715564]
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

popt = [9989.8597311885387, 13254.579666842581, 5179.6582631480323, 16114.516866452037, 4960.964704011526, 10093.046957215538, 1114.7911879070746, 482.75216155154862, 0.00012013273452295402, 1.9328936602124865e-05, 3.5271213084713667e-05, 7.0370582534391241e-05, 4.8628468893893049e-05, 0.00028604144568468854, 9.8687003706937357e-05, 0.00027269742182715564]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 1.326803380978566, 0.51849159072544704, 1.6130874006310818, 0.49660003618702436, 1.0103291966857999, 0.11159227635866342, 0.048324218211431624, 1.2001091670084216, 0.19309336540826341, 0.35235447126515135, 0.70299224872059463, 0.48579158319245663, 2.8575139196964239, 0.98586932432358754, 2.7242089931100368]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9989.85973119,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)