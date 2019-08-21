#current best params = [10059.820072284172, 0.9433998941401327, 15916.518934830679, 23444.504566744272, 0.10749488586606805, 8621.7223991054543, 1846.2952168102534, 1489.2680397615886, 1281.022953784397, 128.75908965247586, 0.0, 0.0, 7.537919926022724e-05, 0.00033859933932101673, 0.00016993327132443421, 0.0, 0.00033436095625619117, 0.0]
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
	growth_func_2 = lambda t: before[1] + (after[1] - before[1]) * (t / T)
	phi = dadi.Integration.two_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, m12=nu33, m21=t2, theta0=theta1)
	before = after
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	before.append((1 - t1) * before[-1])
	before[-2] *= t1]
	T = nu32
	after = [m1_12, m1_21, s1]
	growth_func_1 = lambda t: before[0] + (after[0] - before[0]) * (t / T)
	growth_func_2 = after[1]
	growth_func_3 = lambda t: before[2] * (after[2] / before[2]) ** (t / T)
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.9433998941401327, 1.5821872379887101, 2.330509332998536, 0.10749488586606805, 0.85704538820323206, 0.18353163411908177, 0.14804122032606465, 0.12734054332778233, 0.012799343201696048, 0.0, 0.0, 0.75830118175074224, 3.4062484301637235, 1.7094981338184556, 0.0, 3.3636110591341617, 0.0]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 10059.820072
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
