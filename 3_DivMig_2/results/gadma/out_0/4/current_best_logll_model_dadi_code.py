#current best params = [10282.685925297455, 0.89531913950159459, 13406.805963032688, 24663.347044567236, 0.20913939688360106, 18684.561917344621, 4196.062039365891, 7845.8295852433093, 665.9604305603292, 632.75564502228417, 0.0, 0.0, 4.0891944474860014e-05, 5.2087391685956846e-05, 4.3666771046853256e-05, 0.00052960019796061055, 9.9317789931265349e-05, 0.0006630503318824448]
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
	growth_func_1 = lambda t: before[0] * (after[0] / before[0]) ** (t / T)
	growth_func_2 = lambda t: before[1] * (after[1] / before[1]) ** (t / T)
	growth_func_3 = after[2]
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.89531913950159459, 1.3038233454207986, 2.3985315922069046, 0.20913939688360106, 1.8170896255205924, 0.40807062180541198, 0.76301363692739133, 0.064765221402117698, 0.061536027611772057, 0.0, 0.0, 0.42047902190968811, 0.53559828937464415, 0.44901169204666441, 5.4457125016043166, 1.0212536406578716, 6.8179383154114221]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 10282.685925
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
