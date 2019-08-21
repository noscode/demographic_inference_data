#current best params = [9812.4441703836619, 0.13593999490561492, 32327.677077592605, 7169.326347791417, 0.07867217663034938, 9327.8912333850367, 226954.69498221355, 12090.464918790434, 1724.052437225318, 275.50892419265335, 0.00030622815919306327, 0.0, 0.0, 0.0, 6.7137815969849107e-05, 2.2945988573531492e-06, 0.00014849287969761099, 0.0]
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
	growth_func_2 = after[1]
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

popt = [0.13593999490561492, 3.2945590839808681, 0.73063614154668866, 0.07867217663034938, 0.95061852800537472, 23.129272487197216, 1.2321563016156964, 0.17570061111062693, 0.028077502343830517, 3.0048467154812935, 0.0, 0.0, 0.0, 0.65878607092563701, 0.022515623181203922, 1.4570780917323054, 0.0]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 9812.444170
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
