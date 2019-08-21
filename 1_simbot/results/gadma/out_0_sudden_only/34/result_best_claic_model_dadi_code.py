#current best params = [2271712.2991991378, 2805.2079715639984, 227171229.91991377, 18335.197861319903, 145.16786241047836]
import dadi
import numpy as np

def generated_model((nu21, t1, nu31, t2), ns, pts):
	theta1 = 2.0
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx, theta0=theta1)
	T = nu31
	after = [nu21]
	growth_func = after[0]
	phi = dadi.Integration.one_pop(phi, xx, nu=growth_func,T=T, theta0=theta1)
	T = t2
	after = [t1]
	growth_func = after[0]
	phi = dadi.Integration.one_pop(phi, xx, nu=growth_func,T=T, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*1)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/1_simbot/simulated_data/00.fs')

popt = [0.001234842973977355, 100.0, 0.0080710915144420955, 6.3902397527035174e-05]
ns = [20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 2271712.299199
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
