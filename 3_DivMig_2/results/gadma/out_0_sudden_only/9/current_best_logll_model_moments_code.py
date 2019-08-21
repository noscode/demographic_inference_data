#current best params = [9981.8780651522829, 12832.559705092937, 5190.7349719734311, 16370.792247112891, 4847.0817059653727, 10076.808491022113, 1175.6027835285045, 479.86463531233971, 0.00016878937198129494, 4.2375335409845494e-05, 2.3996015462225361e-05, 5.9884271874819368e-05, 4.5108128120851773e-05, 0.00032735558104791441, 9.4219723219073903e-05, 0.00026260197962113812]
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

popt = [9981.8780651522829, 12832.559705092937, 5190.7349719734311, 16370.792247112891, 4847.0817059653727, 10076.808491022113, 1175.6027835285045, 479.86463531233971, 0.00016878937198129494, 4.2375335409845494e-05, 2.3996015462225361e-05, 5.9884271874819368e-05, 4.5108128120851773e-05, 0.00032735558104791441, 9.4219723219073903e-05, 0.00026260197962113812]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [0.99999999999999989, 1.2855857005399278, 0.52001586656270593, 1.6400513150190579, 0.48558815027875474, 1.0095102770490898, 0.11777370709752999, 0.048073582163620517, 1.6848349298109173, 0.42298543103100761, 0.23952530039344236, 0.59775749987487525, 0.45026383464960923, 3.2676234939673576, 0.94048978850519305, 2.6212609402458056]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9981.87806515,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)