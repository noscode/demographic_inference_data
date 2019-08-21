#current best params = [13086.015364215113, 4804.1309559775173, 6063.5605146492626, 14975.16454332648, 4238.1847960185205, 5364.6790173208219, 18155.136013926953, 1440.2612540709945, 0.0014813596626743193, 0.00038000327065986368, 1.1508475391196158e-05, 0.00012861424333253184, 5.3818303105055383e-05, 0.00044004817899958513, 5.0794012047408137e-05, 0.0014205033275442997]
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

popt = [13086.015364215113, 4804.1309559775173, 6063.5605146492626, 14975.16454332648, 4238.1847960185205, 5364.6790173208219, 18155.136013926953, 1440.2612540709945, 0.0014813596626743193, 0.00038000327065986368, 1.1508475391196158e-05, 0.00012861424333253184, 5.3818303105055383e-05, 0.00044004817899958513, 5.0794012047408137e-05, 0.0014205033275442997]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 0.36711946473139917, 0.46336186729771195, 1.1443639737943008, 0.32387129909752504, 0.40995512140319046, 1.3873693029257637, 0.11006110064714723, 19.385095305684658, 4.97272863830697, 0.15060008578788445, 1.6830479643064127, 0.70426714130874057, 5.7584772313834529, 0.66469122206251041, 18.588728369163398]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=13086.0153642,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)