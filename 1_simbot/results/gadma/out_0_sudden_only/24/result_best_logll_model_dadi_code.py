#current best params = [10042.015898894088, 335.19228707773635, 10453.028136415329, 174.47697661813362, 469.14472807708893]
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

popt = [0.03337898390647346, 1.0409292558047538, 0.017374696313451218, 0.046718182165869218]
ns = [20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 10042.015899
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
