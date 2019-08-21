#current best params = [10131.578956515486, 14014.059714209052, 4086.3216438316913, 15647.757322491079, 4950.6159281128212, 9088.8741661455297, 781.02460109908429, 640.50653317654849, 0.0, 3.0601518251329613e-06, 3.8693022983725177e-05, 7.4331963987926569e-05, 4.3280984317955586e-05, 0.00045060048561770268, 9.0698087724773997e-05, 0.00050557265173540233]
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

popt = [10131.578956515486, 14014.059714209052, 4086.3216438316913, 15647.757322491079, 4950.6159281128212, 9088.8741661455297, 781.02460109908429, 640.50653317654849, 0.0, 3.0601518251329613e-06, 3.8693022983725177e-05, 7.4331963987926569e-05, 4.3280984317955586e-05, 0.00045060048561770268, 9.0698087724773997e-05, 0.00050557265173540233]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 1.3832058926211885, 0.40332525279328069, 1.544453968098251, 0.48863222103492021, 0.89708368312133557, 0.077088142376546112, 0.063218826594116123, 0.0, 0.031004169835259568, 0.39202141742588004, 0.75310016213654374, 0.43850470993307555, 4.5652943978799749, 0.91891483698851573, 5.1222492393121346]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10131.5789565,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)