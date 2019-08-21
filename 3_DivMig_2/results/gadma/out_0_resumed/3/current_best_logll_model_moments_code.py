#current best params = [10104.086393995987, 0.89423067951811241, 8825.0993068863372, 4939.8140232054229, 0.5556815423331023, 15660.323005668113, 4543.7430003799791, 11122.368716921159, 233.382297224578, 1232.8653881878558, 0.0, 0.0, 4.0811544716715676e-05, 5.24102660334341e-05, 3.5504872804179507e-05, 0.00081047460811441954, 9.3697968523441048e-05, 0.00074367689931658157]
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
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: before[1] * (after[1] / before[1]) ** (t / T)]
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

popt = [10104.086393995987, 0.89423067951811241, 8825.0993068863372, 4939.8140232054229, 0.5556815423331023, 15660.323005668113, 4543.7430003799791, 11122.368716921159, 233.382297224578, 1232.8653881878558, 0.0, 0.0, 4.0811544716715676e-05, 5.24102660334341e-05, 3.5504872804179507e-05, 0.00081047460811441954, 9.3697968523441048e-05, 0.00074367689931658157]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.89423067951811241, 0.87341882905220958, 0.48889269455779211, 0.5556815423331023, 1.5498999508727214, 0.44969360149968091, 1.1007792573438655, 0.023097812916886534, 0.12201651293485026, 0.0, 0.0, 0.41236337369012566, 0.52955785593413152, 0.3587443022212683, 8.1891054605281361, 0.94673236890276491, 7.5141756399137956]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10104.086394,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)