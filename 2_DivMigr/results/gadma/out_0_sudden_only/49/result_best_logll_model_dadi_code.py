#current best params = [10004.660638085019, 10004.943173367295, 988.04370052174284, 494.33214500665872, 0.0004993521291752767, 0.00025733989858671723]
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

popt = [1.0000282403664149, 0.098758342363011245, 0.049410186201106196, 4.9958485913038366, 2.5745983539993205]
ns = [20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 10004.660638
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
