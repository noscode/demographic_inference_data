#current best params = [9989.4660170128627, 13260.272953736603, 5179.8271580520213, 16112.333171372218, 4961.5111076850453, 10095.153270706542, 1114.7427994024242, 482.66679968862951, 0.00011989613988540882, 1.9311681448141271e-05, 3.5280234702241428e-05, 7.0423934210191519e-05, 4.8628650245840641e-05, 0.00028585742083221362, 9.8696522105543903e-05, 0.00027235644764340468]
import dadi
import numpy as np

def generated_model((nu21, nu22, t1, m1_12, m1_21, nu31, nu32, nu33, t2, m2_12,
		m2_13, m2_21, m2_23, m2_31, m2_32), ns, pts):
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
	T = nu32
	after = [t1, m1_12, m1_21]
	growth_func_1 = after[0]
	growth_func_2 = after[1]
	growth_func_3 = after[2]
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=m2_12, m13=m2_13, m21=m2_21, m23=m2_23, m31=m2_31, m32=m2_32, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [1.3274256032457885, 0.51852893330137562, 1.6129323773594726, 0.49667430663813195, 1.0105798701866231, 0.11159183058473073, 0.04831757762292891, 1.1976984149563119, 0.192913385557585, 0.35243070563027867, 0.70349749757705782, 0.48577424908402927, 2.855562991114343, 0.98592555357068967, 2.7206954782481341]
ns = [20, 20, 20]
pts = [20, 30, 40]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model = func_ex(popt, ns, pts)
N_A = 9989.466017
ll_model = dadi.Inference.ll(N_A * model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
