#current best params = [10002.463018350449, 0.61084866224050471, 14751.40901937142, 7891.93056672068, 0.20672717184298722, 15547.503264204737, 8632.8754393816289, 15700.996481527502, 1167.5542190262383, 353.68075356284737, 3.4125544697610166e-05, 3.6489682409418377e-05, 4.8111668239359464e-05, 0.00011358293931617212, 4.8042998394693869e-05, 0.00011270040042093606, 9.9578849063998364e-05, 8.2429317159229685e-05]
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
	growth_func_2 = lambda t: before[1] + (after[1] - before[1]) * (t / T)
	growth_func_3 = lambda t: before[2] * (after[2] / before[2]) ** (t / T)
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.61084866224050471, 1.4747776614928332, 0.78899872483829225, 0.20672717184298722, 1.5543674828571119, 0.86307496699001196, 1.5697130249542102, 0.11672667190913392, 0.035359366279484078, 0.34133949881891096, 0.36498669885156021, 0.48123518231533891, 1.1361091500255549, 0.48054831473359538, 1.1272815873637003, 0.99603375517254489, 0.82449619651307504]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 10002.463018
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
