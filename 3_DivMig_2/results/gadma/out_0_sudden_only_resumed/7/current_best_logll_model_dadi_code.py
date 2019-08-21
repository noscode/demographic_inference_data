#current best params = [10022.49, 14707.07007694302, 4446.3390843136176, 15166.561393827767, 5058.0470727455577, 9917.3650173477236, 968.03672998745299, 543.4223675150871, 0.0, 8.9648302589830816e-05, 4.7279542926963532e-05, 9.0618210262992068e-05, 3.4995573579229827e-05, 0.00035187383071244279, 8.129780817389596e-05, 0.00037068468589426557]
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

popt = [1.4674068097791089, 0.44363617068349459, 1.5132528337596514, 0.50466970510776843, 0.98951109129046011, 0.096586450072532179, 0.05422029530736245, 0.0, 0.89849921622355344, 0.47385874619006274, 0.90822010617873539, 0.35074278624209515, 3.5266519495771509, 0.81480646944479052, 3.7151835575284178]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 10022.490000
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
