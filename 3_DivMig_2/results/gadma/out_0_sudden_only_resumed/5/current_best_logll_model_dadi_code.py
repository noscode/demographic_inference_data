#current best params = [10053.244189275547, 14453.583806269016, 4604.8713085949121, 15431.194685347311, 4931.5251413773713, 9608.8167927083941, 900.87406195502808, 562.62477984261125, 4.6349800328939502e-06, 0.0, 4.5947631337411954e-05, 8.245168724271601e-05, 4.6184832359399291e-05, 0.00040478962857754474, 9.5402759699420226e-05, 0.00039548855646294554]
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

popt = [1.4377034451911155, 0.4580482898751459, 1.5349467689055813, 0.49054067010907298, 0.95579263885370924, 0.089610283505900445, 0.055964499543619949, 0.046596586083099287, 0.0, 0.46192275775381175, 0.82890694566879952, 0.46430739754979616, 4.0694489813762083, 0.95910723958904776, 3.9759430321860814]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 10053.244189
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
