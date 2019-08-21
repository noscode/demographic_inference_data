#current best params = [9974.0779455657648, 0.52025567263296191, 32961.201033071236, 6434.9641708275303, 0.43108558695524141, 4011.8266001295565, 6637.6341391669057, 16004.498364563524, 1108.3437908154863, 413.82239800641577, 0.0, 0.0, 6.5461419018799966e-05, 0.00013331232663997249, 5.6538360940587944e-05, 0.00020838267700006756, 0.00011547207046424264, 0.00025633596876825789]
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
	growth_func_2 = lambda t: before[1] * (after[1] / before[1]) ** (t / T)
	phi = dadi.Integration.two_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, m12=nu33, m21=t2, theta0=theta1)
	before = after
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	before.append((1 - t1) * before[-1])
	before[-2] *= t1]
	T = nu32
	after = [m1_12, m1_21, s1]
	growth_func_1 = lambda t: before[0] + (after[0] - before[0]) * (t / T)
	growth_func_2 = lambda t: before[1] + (after[1] - before[1]) * (t / T)
	growth_func_3 = lambda t: before[2] + (after[2] - before[2]) * (t / T)
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.52025567263296191, 3.3046865297182677, 0.64516882722862223, 0.43108558695524141, 0.40222531065271233, 0.66548849682068489, 1.604609313453254, 0.11112243125272842, 0.041489789859762552, 0.0, 0.0, 0.652917295720852, 1.329667537011809, 0.56391801893595506, 2.0784250629043282, 1.1517274313462185, 2.5567149327467158]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 9974.077946
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
