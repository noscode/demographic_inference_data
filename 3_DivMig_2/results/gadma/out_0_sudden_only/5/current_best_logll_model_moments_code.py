#current best params = [10085.16899949125, 14359.18867587338, 4342.8507868920933, 15493.94243673032, 4854.1495940897312, 9468.9943510997109, 839.37256235680331, 599.21546019635207, 0.0, 0.0, 4.3725852209151167e-05, 7.6352556019732674e-05, 4.4292516725411491e-05, 0.00046664137034834526, 9.3458395170886016e-05, 0.00042362091550696672]
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

popt = [10085.16899949125, 14359.18867587338, 4342.8507868920933, 15493.94243673032, 4854.1495940897312, 9468.9943510997109, 839.37256235680331, 599.21546019635207, 0.0, 0.0, 4.3725852209151167e-05, 7.6352556019732674e-05, 4.4292516725411491e-05, 0.00046664137034834526, 9.3458395170886016e-05, 0.00042362091550696672]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 1.423792568731147, 0.43061755208179159, 1.5363096481092104, 0.48131564224006562, 0.93890289310743114, 0.083228408210030566, 0.059415510064985494, 0.0, 0.0, 0.44098260917606735, 0.77002843100212703, 0.44669751658856766, 4.7061570821172474, 0.94254370971962242, 4.2722885246069628]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10085.1689995,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)