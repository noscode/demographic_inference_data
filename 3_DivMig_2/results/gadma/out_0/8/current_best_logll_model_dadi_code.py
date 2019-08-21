#current best params = [10086.398888810971, 0.81506487452608045, 25869.851089414617, 11576.762120521949, 0.14134939066480093, 12211.277376966575, 10305.260416915738, 8965.0638321024835, 1137.4641084332629, 297.42207877457304, 7.0461730667803908e-06, 2.2379855981683775e-05, 6.0570166862216913e-05, 0.00014706694775333643, 5.7764428900004648e-05, 0.00016559145889449693, 0.00012177251199631997, 8.5535361097025917e-06]
import dadi
import numpy as np

def generated_model((s0, nu21, nu22, t1, m1_12, m1_21, s1, nu31, nu32, nu33,
		t2, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32), ns, pts):
	theta1 = 2.0
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx, theta0=theta1)
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	before = [s0, 1 - s0]
	T = nu31
	after = [nu21, nu22]
	growth_func_1 = lambda t: before[0] + (after[0] - before[0]) * (t / T)
	growth_func_2 = lambda t: before[1] + (after[1] - before[1]) * (t / T)
	phi = dadi.Integration.two_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, m12=nu33, m21=t2, theta0=theta1)
	before = after
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	before.append((1 - t1) * before[-1])
	before[-2] *= t1]
	T = nu32
	after = [m1_12, m1_21, s1]
	growth_func_1 = after[0]
	growth_func_2 = lambda t: before[1] * (after[1] / before[1]) ** (t / T)
	growth_func_3 = after[2]
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.81506487452608045, 2.5648253033213391, 1.1477596958181244, 0.14134939066480093, 1.2106677032684847, 1.0216986786381763, 0.88882701655271579, 0.11277207266659589, 0.029487439675274875, 0.071070512191143534, 0.22573215450540479, 0.61093486373415973, 1.4833758984000738, 0.5826350714698072, 1.6702215069900415, 1.2282461296874025, 0.086274377112308745]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 10086.398889
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
