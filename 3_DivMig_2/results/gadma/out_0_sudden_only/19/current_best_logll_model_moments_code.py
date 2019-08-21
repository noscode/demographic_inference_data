#current best params = [9966.5266487159843, 12978.326261081078, 5348.5381388481173, 15980.712278915997, 4957.6719563598626, 9921.9457897074517, 1178.7771760806104, 473.82362424774078, 0.00016027807123787896, 1.4782847613878896e-05, 3.1027064033337124e-05, 6.2509173140633972e-05, 5.0420882821972815e-05, 0.00027415803328018523, 0.00010332733015642413, 0.00026989185296157233]
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

popt = [9966.5266487159843, 12978.326261081078, 5348.5381388481173, 15980.712278915997, 4957.6719563598626, 9921.9457897074517, 1178.7771760806104, 473.82362424774078, 0.00016027807123787896, 1.4782847613878896e-05, 3.1027064033337124e-05, 6.2509173140633972e-05, 5.0420882821972815e-05, 0.00027415803328018523, 0.00010332733015642413, 0.00026989185296157233]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 1.3021914974514328, 0.53665016182314473, 1.6034384738212522, 0.49743226814113556, 0.99552694127253694, 0.11827361904786314, 0.047541499756967469, 1.5974156681971197, 0.14733364468763152, 0.30923206051967173, 0.62299933989532996, 0.50252107229697807, 2.7324033446465297, 1.0298145895446758, 2.6898843448128464]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9966.52664872,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)