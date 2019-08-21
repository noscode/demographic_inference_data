#current best params = [9939.9933437368491, 0.50585529227157744, 14407.492410550451, 6995.4615283640651, 0.29200374468281681, 15350.426577081484, 10935.296903735825, 22331.236852558319, 1239.3421350874137, 342.90541171478691, 7.2340883890521256e-05, 2.056148887051433e-05, 4.1617481772751675e-05, 0.00010177908615700289, 5.2407861435677078e-05, 7.7476100555112009e-05, 0.00010678634892277276, 1.7499066760357038e-05]
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
	growth_func_1 = after[0]
	growth_func_2 = lambda t: before[1] * (after[1] / before[1]) ** (t / T)
	phi = dadi.Integration.two_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, m12=nu33, m21=t2, theta0=theta1)
	before = after
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	before.append((1 - t1) * before[-1])
	before[-2] *= t1]
	T = nu32
	after = [m1_12, m1_21, s1]
	growth_func_1 = after[0]
	growth_func_2 = lambda t: before[1] * (after[1] / before[1]) ** (t / T)
	growth_func_3 = lambda t: before[2] * (after[2] / before[2]) ** (t / T)
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.50585529227157744, 1.4494468871680437, 0.70376923670395397, 0.29200374468281681, 1.5443095428986104, 1.1001312099093219, 2.2466048095121858, 0.12468239084569592, 0.034497549430538628, 0.71906790435182155, 0.20438106251023175, 0.41367749180424129, 1.011683438932228, 0.52093379383011329, 0.77011192381650018, 1.0614555974943218, 0.1739406071195557]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 9939.993344
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
