import matplotlib
matplotlib.use('Agg')
import moments
def model_1((nuB, nuF, tB, tF), ns):
        sts = moments.LinearSystem_1D.steady_state_1D(ns[0])
        fs = moments.Spectrum(sts)
        fs.integrate([nuB], tB)
        fs.integrate([nuF], tF)
        return fs

theta = 20000
popt = [0.01, 1.0, 0.005, 0.05]
ns  = [20]

data = model_1(popt, ns) * theta
model = data
ll_model = moments.Inference.ll_multinom(model, model)
print('Maximum log composite likelihood: {0}'.format(ll_model))

model = moments.ModelPlot.generate_model(model_1, popt, ns)
moments.ModelPlot.plot_model(model, 
	save_file='model_1_pop.png',
	fig_title='Maximum log composite likelihood: %.2f' % ll_model,
	pop_labels=['Pop 1'],
	nref=10000,
	draw_scale=False,
	gen_time=1.0,
	gen_time_units="Genetic units",
	reverse_timeline=True)
