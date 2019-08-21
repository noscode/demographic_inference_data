#current best params = [10038.575526696555, 14647.557416348116, 4783.712374548898, 0.46413502239506038, 15851.153465447318, 4936.8618952510897, 9861.2239557485827, 948.31360760844404, 531.73294645806902, 0.0, 0.0, 4.9978278105055886e-05, 9.0331951450151019e-05, 4.7438562053327608e-05, 0.00037044936128920576, 9.7366140484361898e-05, 0.00033842463202174992]
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
	growth_func_1 = lambda t: before[0] * (after[0] / before[0]) ** (t / T)
	growth_func_2 = after[1]
	growth_func_3 = after[2]
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [1.4591270820640287, 0.47653298636117331, 0.46413502239506038, 1.5790241776129303, 0.4917890872188006, 0.98233299431016641, 0.094466949527500385, 0.052968964077022694, 0.0, 0.0, 0.50171071945184831, 0.90680411710622744, 0.47621558805021041, 3.7187838921181915, 0.97741735499521409, 3.3973012286448263]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 10038.575527
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
