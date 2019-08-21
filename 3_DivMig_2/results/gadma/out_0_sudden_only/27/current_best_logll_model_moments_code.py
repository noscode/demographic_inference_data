#current best params = [9865.5430998179418, 14771.472857649289, 6077.9506804859848, 14402.871115064147, 4683.9754305188972, 10497.401892584081, 1274.4326330444753, 354.447309396582, 7.6645404744126217e-05, 0.0, 9.9096618895012862e-12, 0.00015781671548325876, 6.7733634109315915e-05, 4.3722245599359571e-06, 0.00011560791994617855, 0.0]
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

popt = [9865.5430998179418, 14771.472857649289, 6077.9506804859848, 14402.871115064147, 4683.9754305188972, 10497.401892584081, 1274.4326330444753, 354.447309396582, 7.6645404744126217e-05, 0.0, 9.9096618895012862e-12, 0.00015781671548325876, 6.7733634109315915e-05, 4.3722245599359571e-06, 0.00011560791994617855, 0.0]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 1.497279238273449, 0.61607867088413482, 1.4599166988921204, 0.47478130530952067, 1.064047035867473, 0.12918017996069509, 0.035927805069659362, 0.7561485439061677, 0.0, 9.7764196475498235e-08, 1.5569476084717948, 0.66822908661275482, 0.043134369838130716, 1.1405349169093268, 0.0]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9865.54309982,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)