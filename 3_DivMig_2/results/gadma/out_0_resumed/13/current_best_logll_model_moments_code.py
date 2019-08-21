#current best params = [10004.884123469577, 0.12543216486423883, 19813.917129352943, 3334.2999807370602, 0.57972375655833963, 13616.278339394097, 4738.9925674570304, 11643.55163968436, 246.56625229672946, 1371.3923429781776, 0.00018903734878048662, 0.00020029459159407549, 4.8701771299017758e-05, 6.0888272830147306e-05, 3.8967493862144323e-05, 0.00079031037844855723, 9.727953123604356e-05, 0.00080210451906456942]
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
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: after[1]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, Ms[0]],[Ms[1], 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - Ns[4]) * before[-1])
	before[-2] *= Ns[4]
	T = Ts[1]
	after = Ns[5:8]
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T), lambda t: before[1] + (after[1] - before[1]) * (t / T), lambda t: before[2] * (after[2] / before[2]) ** (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [10004.884123469577, 0.12543216486423883, 19813.917129352943, 3334.2999807370602, 0.57972375655833963, 13616.278339394097, 4738.9925674570304, 11643.55163968436, 246.56625229672946, 1371.3923429781776, 0.00018903734878048662, 0.00020029459159407549, 4.8701771299017758e-05, 6.0888272830147306e-05, 3.8967493862144323e-05, 0.00079031037844855723, 9.727953123604356e-05, 0.00080210451906456942]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [0.99999999999999989, 0.12543216486423883, 1.9804244491820966, 0.33326722624557131, 0.57972375655833963, 1.360963122746506, 0.47366791148937393, 1.1637867561475075, 0.024644588508359769, 0.13707228650066508, 1.8912967695566718, 2.0039241794563889, 0.4872555784543891, 0.60918011414382489, 0.38986526067276595, 7.9069637579532044, 0.97327043760205512, 8.0249627681523119]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10004.8841235,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)