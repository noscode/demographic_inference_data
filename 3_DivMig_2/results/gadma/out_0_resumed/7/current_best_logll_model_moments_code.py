#current best params = [16437.57576974465, 0.26185605111073817, 201.96913916934525, 10168.031482352897, 0.03693520965488966, 17323.020809559282, 3193.496889172457, 5309.5946340396094, 21269.597348858755, 11102.612139861732, 0.00020357181484132169, 6.8901091493163223e-06, 0.00027385160540629571, 4.0768157965384667e-07, 0.0, 0.00139651274355663, 0.00016313114930937859, 0.00077561089090817442]
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
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T), lambda t: before[1] * (after[1] / before[1]) ** (t / T), lambda t: after[2]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [16437.57576974465, 0.26185605111073817, 201.96913916934525, 10168.031482352897, 0.03693520965488966, 17323.020809559282, 3193.496889172457, 5309.5946340396094, 21269.597348858755, 11102.612139861732, 0.00020357181484132169, 6.8901091493163223e-06, 0.00027385160540629571, 4.0768157965384667e-07, 0.0, 0.00139651274355663, 0.00016313114930937859, 0.00077561089090817442]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.26185605111073817, 0.012287039281126474, 0.61858461520027741, 0.03693520965488966, 1.0538671305439336, 0.1942802840215937, 0.32301567508589446, 1.2939619349471245, 0.67544097106444589, 3.3462271310386535, 0.1132566912036979, 4.5014565135321991, 0.0067012968554892929, 0.0, 22.955284035626086, 2.6814806271784382, 12.749162787142268]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=16437.5757697,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)