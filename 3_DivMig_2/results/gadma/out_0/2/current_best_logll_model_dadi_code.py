#current best params = [10036.254121378481, 0.61078545220971847, 23557.967544906998, 6597.3664270557028, 0.0015212515830892686, 9977.3334633221493, 4801.4255994170062, 12061.925351736858, 1126.7733259928282, 445.80430102291018, 0.0001163438400961223, 2.8207902332132722e-05, 3.428474789101832e-05, 7.2351517048459277e-05, 4.9071597282860627e-05, 0.00024839664866595538, 0.00010449686810685563, 0.0002913370790100933]
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
	growth_func_1 = lambda t: before[0] * (after[0] / before[0]) ** (t / T)
	growth_func_2 = after[1]
	growth_func_3 = lambda t: before[2] + (after[2] - before[2]) * (t / T)
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.61078545220971847, 2.3472868721733113, 0.65735346547298801, 0.0015212515830892686, 0.99412921819797062, 0.47840813328843151, 1.2018353865754996, 0.11227030646749564, 0.044419391501186789, 1.1676563446617063, 0.2831016760363087, 0.34409044232165481, 0.72613821116558486, 0.49249502047273502, 2.4929718889102972, 1.0487571232085733, 2.9239329599254167]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 10036.254121
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
