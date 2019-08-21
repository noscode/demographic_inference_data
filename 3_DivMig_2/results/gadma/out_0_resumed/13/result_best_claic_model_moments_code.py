#current best params = [10004.884123469576, 0.12543216486423883, 19813.91712935294, 3334.2999807370602, 0.57972375655833963, 13616.278339394095, 4738.9925674570304, 11643.551639684358, 246.56625229672946, 1371.3923429781773, 0.00018903734878048664, 0.00020029459159407549, 4.8701771299017758e-05, 6.0888272830147313e-05, 3.8967493862144323e-05, 0.00079031037844855734, 9.7279531236043573e-05, 0.00080210451906456942]
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

popt = [10004.884123469576, 0.12543216486423883, 19813.91712935294, 3334.2999807370602, 0.57972375655833963, 13616.278339394095, 4738.9925674570304, 11643.551639684358, 246.56625229672946, 1371.3923429781773, 0.00018903734878048664, 0.00020029459159407549, 4.8701771299017758e-05, 6.0888272830147313e-05, 3.8967493862144323e-05, 0.00079031037844855734, 9.7279531236043573e-05, 0.00080210451906456942]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.12543216486423883, 1.9804244491820968, 0.33326722624557137, 0.57972375655833963, 1.3609631227465062, 0.47366791148937404, 1.1637867561475077, 0.024644588508359776, 0.13707228650066511, 1.8912967695566716, 2.0039241794563885, 0.48725557845438899, 0.60918011414382478, 0.38986526067276589, 7.9069637579532026, 0.97327043760205501, 8.0249627681523101]
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