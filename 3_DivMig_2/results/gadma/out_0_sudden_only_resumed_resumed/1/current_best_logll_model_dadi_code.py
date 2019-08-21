#current best params = [9989.47, 15082.86, 4765.51, 14862.84, 5075.75, 9988.48, 1012.54, 512.47, 4.65e-06, 5.2e-05, 5.22e-05, 9.44e-05, 4.11e-05, 0.0003, 8.98e-05, 0.000328]
import dadi
import numpy as np

def generated_model((nu21, nu22, t1, m1_12, m1_21, nu31, nu32, nu33, t2, m2_12,
		m2_13, m2_21, m2_23, m2_31, m2_32), ns, pts):
	theta1 = 2.0
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx, theta0=theta1)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	T = nu31
	after = [nu21, nu22]
	growth_func_1 = after[0]
	growth_func_2 = after[1]
	phi = dadi.Integration.two_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, m12=nu33, m21=t2, theta0=theta1)
	before = after
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	T = nu32
	after = [t1, m1_12, m1_21]
	growth_func_1 = after[0]
	growth_func_2 = after[1]
	growth_func_3 = after[2]
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [1.509875899321986, 0.47705333716403375, 1.4878507067942544, 0.5081100398719852, 0.9999008956431122, 0.10136073285169284, 0.05130101997403266, 0.0464510355, 0.5194524399999999, 0.521450334, 0.943005968, 0.410567217, 2.9968409999999994, 0.897054406, 3.27654616]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 9989.470000
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
