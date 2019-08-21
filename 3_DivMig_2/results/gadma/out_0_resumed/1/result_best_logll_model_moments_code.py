#current best params = [10004.331677259052, 0.5602902407505379, 23036.625473923228, 8318.1923527994168, 0.21486758723771657, 7801.6531106646307, 11006.230543470137, 12868.300464671882, 1319.3079626276512, 323.29397420585622, 0.00019695746443910262, 1.8439921183517536e-06, 2.9052664545407898e-05, 5.9799340241315932e-06, 6.5881667535556265e-05, 0.00012683838335268157, 0.00013268291253953934, 5.0433792685190809e-06]
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
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: before[1] * (after[1] / before[1]) ** (t / T), lambda t: before[2] + (after[2] - before[2]) * (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [10004.331677259052, 0.5602902407505379, 23036.625473923228, 8318.1923527994168, 0.21486758723771657, 7801.6531106646307, 11006.230543470137, 12868.300464671882, 1319.3079626276512, 323.29397420585622, 0.00019695746443910262, 1.8439921183517536e-06, 2.9052664545407898e-05, 5.9799340241315932e-06, 6.5881667535556265e-05, 0.00012683838335268157, 0.00013268291253953934, 5.0433792685190809e-06]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.5602902407505379, 2.3026651071843225, 0.83145907404365504, 0.21486758723771657, 0.77982751495521163, 1.1001465063866798, 1.2862728745712166, 0.13187367284379262, 0.03231539943250173, 1.9704278005607374, 0.01844790876224247, 0.29065249222040518, 0.059825243385538893, 0.65910205347661477, 1.2689332564675593, 1.3274038649503055, 0.050455638976477021]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10004.3316773,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)