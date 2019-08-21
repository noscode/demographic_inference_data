#current best params = [9988.7578109065798, 0.99588001801329151, 22922.763347972621, 5453.6965795663864, 0.51112763450575449, 10060.084071640145, 6676.0831241868955, 17512.874390622146, 1064.6157792118693, 459.00999369789952, 1.7374757628513582e-05, 0.0, 5.5080736745953336e-05, 0.00010386915919549818, 5.0410297780253953e-05, 0.000274722818052805, 0.00010508712803047533, 0.00033754041901387595]
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
	growth_func_1 = lambda t: before[0] * (after[0] / before[0]) ** (t / T)
	growth_func_2 = after[1]
	phi = dadi.Integration.two_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, m12=nu33, m21=t2, theta0=theta1)
	before = after
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	before.append((1 - t1) * before[-1])
	before[-2] *= t1]
	T = nu32
	after = [m1_12, m1_21, s1]
	growth_func_1 = lambda t: before[0] * (after[0] / before[0]) ** (t / T)
	growth_func_2 = lambda t: before[1] + (after[1] - before[1]) * (t / T)
	growth_func_3 = lambda t: before[2] + (after[2] - before[2]) * (t / T)
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [0.99588001801329151, 2.2948562555940226, 0.54598346288980737, 0.51112763450575449, 1.0071406537313063, 0.66835969502607995, 1.7532584854044706, 0.10658139874504022, 0.045952660219343106, 0.17355224597442373, 0.0, 0.55018813940163047, 1.0375238752263314, 0.50353625570263827, 2.7441396946592231, 1.0496898709401503, 3.3716094969215331]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 9988.757811
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
