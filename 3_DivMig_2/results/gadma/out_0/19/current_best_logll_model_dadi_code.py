#current best params = [9967.3146445165075, 15330.245064992352, 5041.7803921540635, 0.80046644362610975, 14072.792665630384, 5768.805130313388, 10599.895470591055, 1070.4862053025111, 462.62585227749247, 2.0353352299054581e-10, 4.5591126143998023e-05, 4.7862183191897412e-05, 0.00011366669235346279, 4.388469326274478e-05, 0.00026237246902314663, 8.9994004175019657e-05, 0.00020443544383591698]
import dadi
import numpy as np

def generated_model((nu21, nu22, t1, m1_12, m1_21, s1, nu31, nu32, nu33, t2,
		m2_12, m2_13, m2_21, m2_23, m2_31, m2_32), ns, pts):
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
	before.append((1 - t1) * before[-1])
	before[-2] *= t1]
	T = nu32
	after = [m1_12, m1_21, s1]
	growth_func_1 = lambda t: before[0] + (after[0] - before[0]) * (t / T)
	growth_func_2 = lambda t: before[1] + (after[1] - before[1]) * (t / T)
	growth_func_3 = after[2]
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [1.5380516831006481, 0.50583136701998122, 0.80046644362610975, 1.4118940925951902, 0.57877225070717331, 1.0634655219219511, 0.10739966013729046, 0.046414291991073532, 2.0286826643537044e-06, 0.45442109927507091, 0.47705743944713092, 1.1329516872884222, 0.43741254582787098, 2.6151489528123633, 0.8969985557323531, 2.0376723932039673]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 9967.314645
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
