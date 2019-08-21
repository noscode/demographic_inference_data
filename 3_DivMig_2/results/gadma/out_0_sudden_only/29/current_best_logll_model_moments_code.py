#current best params = [9845.1209359532404, 93851.252613459263, 4184.8284910254724, 5900.1726280049143, 4952.0078513346307, 11290.438575056833, 1221.4306183828669, 366.31524901795143, 0.0, 0.00028457890089179846, 4.5375404451365484e-05, 0.00011441288541497902, 2.8185550585132386e-05, 0.0, 7.5210908736372508e-05, 2.0696485350407467e-05]
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

popt = [9845.1209359532404, 93851.252613459263, 4184.8284910254724, 5900.1726280049143, 4952.0078513346307, 11290.438575056833, 1221.4306183828669, 366.31524901795143, 0.0, 0.00028457890089179846, 4.5375404451365484e-05, 0.00011441288541497902, 2.8185550585132386e-05, 0.0, 7.5210908736372508e-05, 2.0696485350407467e-05]
ns = [20, 20, 20]
model = generated_model(popt, ns)
ll_model = moments.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
#now we need to norm vector of params so that N_A is 1:
popt_norm = [1.0, 9.5327678780181735, 0.4250662351686269, 0.59929915197467687, 0.5029910636496574, 1.1468054733411612, 0.12406456216523902, 0.037207795760051114, 0.0, 2.8017136951004074, 0.44672634434148417, 1.1264086935418292, 0.27749015415705597, 0.0, 0.74046049221252941, 0.20375940122394609]
print('Drawing model to model_from_GADMA.png')
model = moments.ModelPlot.generate_model(generated_model, popt_norm, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_from_GADMA.png',
	fig_title='Demographic model from GADMA',
	pop_labels=['pop 0', 'pop 1', 'pop 2'],
	nref=9845.12093595,
	draw_scale=True,
	gen_time=1.0,
	gen_time_units='Genetic units',
	reverse_timeline=True)