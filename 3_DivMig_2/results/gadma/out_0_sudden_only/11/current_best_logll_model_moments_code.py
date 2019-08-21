#current best params = [9969.0796549828519, 12527.834749857138, 5425.9090384443225, 16724.766261291315, 4944.5684442753, 10133.217691002747, 1217.8013607374576, 462.5621747309529, 0.00019806140219660859, 1.9876696690564564e-05, 2.6075155282141072e-05, 4.487405024164919e-05, 5.0310599952862643e-05, 0.00024942830637516151, 0.00010124053339284976, 0.00024646028578512668]
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

popt = [9969.0796549828519, 12527.834749857138, 5425.9090384443225, 16724.766261291315, 4944.5684442753, 10133.217691002747, 1217.8013607374576, 462.5621747309529, 0.00019806140219660859, 1.9876696690564564e-05, 2.6075155282141072e-05, 4.487405024164919e-05, 5.0310599952862643e-05, 0.00024942830637516151, 0.00010124053339284976, 0.00024646028578512668]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 1.2566691393217368, 0.54427381726579815, 1.6776640211648588, 0.49599046405490932, 1.0164647130628406, 0.12215785236792276, 0.046399686905877029, 1.9744898950755867, 0.19815237258617219, 0.25994530002371119, 0.44735298130070328, 0.50155037842006422, 2.4865706544614521, 1.0092749417062705, 2.4569822207817658]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9969.07965498,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)