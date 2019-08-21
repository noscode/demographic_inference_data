#current best params = [10259.599186780875, 0.94053997071489648, 13531.544544266328, 2799.1761804826674, 0.9406773714537684, 16989.422266715221, 4398.7033730281573, 7393.1388951379067, 223.55887949488533, 1088.5534599168054, 0.0, 0.0, 3.5880514938518348e-05, 3.2866290261737463e-05, 3.6086456254305522e-05, 0.0007289254155500133, 8.7517267306153933e-05, 0.00075352287794549328]
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
	growth_funcs = [lambda t: before[0] * (after[0] / before[0]) ** (t / T), lambda t: before[1] + (after[1] - before[1]) * (t / T)]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0, Ms[0]],[Ms[1], 0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	before = after
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
	before.append((1 - Ns[4]) * before[-1])
	before[-2] *= Ns[4]
	T = Ts[1]
	after = Ns[5:8]
	growth_funcs = [lambda t: before[0] + (after[0] - before[0]) * (t / T), lambda t: after[1], lambda t: after[2]]
	list_growth_funcs = lambda t: [ f(t) for f in growth_funcs]
	m = np.array([[0.0, Ms[2], Ms[3]],[Ms[4], 0.0, Ms[5]], [Ms[6], Ms[7], 0.0]])
	fs.integrate(Npop=list_growth_funcs, tf=T, m=m, dt_fac=0.1, theta=theta1)

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [10259.599186780875, 0.94053997071489648, 13531.544544266328, 2799.1761804826674, 0.9406773714537684, 16989.422266715221, 4398.7033730281573, 7393.1388951379067, 223.55887949488533, 1088.5534599168054, 0.0, 0.0, 3.5880514938518348e-05, 3.2866290261737463e-05, 3.6086456254305522e-05, 0.0007289254155500133, 8.7517267306153933e-05, 0.00075352287794549328]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.94053997071489648, 1.3189155149161418, 0.27283484759222421, 0.9406773714537684, 1.6559538006714221, 0.4287402746391622, 0.72060699063796763, 0.021790215721383434, 0.10610097335180184, 0.0, 0.0, 0.36811970188450188, 0.33719496484182587, 0.37023257724047653, 7.4784826006008283, 0.89789208448350133, 7.7308427057903675]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=10259.5991868,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)