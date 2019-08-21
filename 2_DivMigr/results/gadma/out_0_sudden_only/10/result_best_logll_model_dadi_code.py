#current best params = [9997.8253298638519, 9982.717510400953, 1012.1371499429113, 505.58986220152758, 0.00050242145472506735, 0.00024486964410005789]
import dadi
import numpy as np

def generated_model((nu21, nu22, t1, m1_12, m1_21), ns, pts):
	theta1 = 2.0
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx, theta0=theta1)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	T = t1
	after = [nu21, nu22]
	growth_func_1 = after[0]
	growth_func_2 = after[1]
	phi = dadi.Integration.two_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, m12=m1_12, m21=m1_21, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*2)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/2_DivMigr/simulated_data/00.fs')

popt = [0.99848888943700875, 0.10123573042625804, 0.050569983523448152, 5.0231219463173229, 2.4481639302983051]
ns = [20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 9997.825330
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
