#current best params = [9999.6952524595272, 9989.2163713949685, 1004.9638000858441, 503.30917607321078, 0.00050125585730639367, 0.00025038083056456957]
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

popt = [0.9989520799584386, 0.10049944270438271, 0.050332451476400419, 5.0124058165742751, 2.5037320027033996]
ns = [20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 9999.695252
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
