#current best params = [9999.7612691569229, 9966.6843743269674, 1011.7974064101976, 507.71098941664388, 0.00050275461628207974, 0.00024870423410106658]
import matplotlib
matplotlib.use("Agg")
import moments
import numpy as np

def generated_model(params, ns):
	Ns = params[:3]
	Ts = params[3:4]
	Ms = params[4:]
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

	return fs
data = moments.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/2_DivMigr/simulated_data/00.fs')

popt = [9999.7612691569229, 9966.6843743269674, 1011.7974064101976, 507.71098941664388, 0.00050275461628207974, 0.00024870423410106658]
ns = [20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.99669223155036946, 0.10118215617116447, 0.050772311033326183, 5.0274261397873916, 2.4869829676391819]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1'],
	nref=9999.76126916,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)