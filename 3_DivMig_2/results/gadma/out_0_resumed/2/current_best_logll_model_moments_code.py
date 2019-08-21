#current best params = [10506.546802449469, 0.86311348548382782, 7268.6986448789175, 4267.2619334231522, 0.63973863332418246, 17232.643702446803, 4329.5184182676594, 8804.6810798981169, 6586.1526282883542, 1645.4973351096537, 0.00069343353622356981, 0.00094144257205777682, 0.0, 0.00019467918244123203, 5.0158671637143555e-05, 0.0007904156202051645, 8.5427546839473301e-05, 0.00094144257205777682]
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
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T), lambda t: before[1] * (after[1] / before[1]) ** (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, Ms[0]],[Ms[1], 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - Ns[4]) * before[-1])
	before[-2] *= Ns[4]
	T = Ts[1]
	after = Ns[5:8]
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: before[1] + (after[1] - before[1]) * (t / T), lambda t: before[2] * (after[2] / before[2]) ** (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [10506.546802449469, 0.86311348548382782, 7268.6986448789175, 4267.2619334231522, 0.63973863332418246, 17232.643702446803, 4329.5184182676594, 8804.6810798981169, 6586.1526282883542, 1645.4973351096537, 0.00069343353622356981, 0.00094144257205777682, 0.0, 0.00019467918244123203, 5.0158671637143555e-05, 0.0007904156202051645, 8.5427546839473301e-05, 0.00094144257205777682]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.86311348548382782, 0.69182565704502563, 0.4061526602088037, 0.63973863332418246, 1.6401815007790408, 0.41207815466622072, 0.83801854647860285, 0.62686177981455105, 0.15661638081943599, 7.2855919027209755, 9.8913104451434393, 0.0, 2.0454059417814032, 0.52699443110434352, 8.3045387070726857, 0.89754851908737043, 9.8913104451434393]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10506.5468024,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)