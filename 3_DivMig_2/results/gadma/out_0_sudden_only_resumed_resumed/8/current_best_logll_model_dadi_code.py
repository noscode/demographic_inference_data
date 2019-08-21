#current best params = [9987.750967291764, 11684.893446655095, 5377.9144526793198, 17668.440531745444, 4962.8662964980858, 9799.886461551072, 1249.2959503263428, 475.50714832073965, 0.00024473717513494558, 2.8458831013344049e-05, 2.2943790522901124e-05, 2.9789623547302102e-05, 4.8075673021211369e-05, 0.00025795700223100283, 0.0001000078995202745, 0.00030231006756054166]
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

popt = [1.1699223864232491, 0.53845099565368637, 1.7690109204371105, 0.49689527830145686, 0.98119050961964138, 0.12508280937496147, 0.047609031290222099, 2.4443739576863064, 0.2842397169815199, 0.22915686598844531, 0.29753134139982407, 0.48016784972080639, 2.576410298552382, 0.99885399517043916, 3.0193976696998384]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 9987.750967
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
