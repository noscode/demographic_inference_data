#current best params = [9996.9882692370065, 0.48667262701538139, 14693.7077739742, 4661.0841320895515, 0.54289634907552176, 15481.018211297387, 4972.2187531479649, 10142.326190192493, 1012.621416177032, 512.63940465703138, 2.3143394563933557e-05, 2.7811199599886778e-05, 4.5270475247589694e-05, 9.3721328691165031e-05, 4.4757048960255865e-05, 0.00034067378168311287, 9.3019688674096033e-05, 0.00029748422713635892]
import dadi
import numpy as np

def generated_model(params, ns, pts):
	Ns = params[:8]
	Ts = params[8:10]
	Ms = params[10:]
	theta1 = 2.0
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx, theta0=theta1, nu=Ns[0])
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	before.append((1 - Ns[1]) * before[-1])
	before[-2] *= Ns[1]
	T = Ts[0]
	after = Ns[2:4]
	growth_func_1 = after[0]
	growth_func_2 = lambda t: before[1] * (after[1] / before[1]) ** (t / T)
	phi = dadi.Integration.two_pops(phi, xx,  T=T,nu1=growth_func_1, nu2=growth_func_2, m12=Ms[0], m21=Ms[1], theta0=theta1)
	before = after
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	before.append((1 - Ns[4]) * before[-1])
	before[-2] *= Ns[4]
	T = Ts[1]
	after = Ns[5:8]
	growth_func_1 = lambda t: before[0] + (after[0] - before[0]) * (t / T)
	growth_func_2 = after[1]
	growth_func_3 = after[2]
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=Ms[2], m13=Ms[3], m21=Ms[4], m23=Ms[5], m31=Ms[6], m32=Ms[7], theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [9996.9882692370065, 0.48667262701538139, 14693.7077739742, 4661.0841320895515, 0.54289634907552176, 15481.018211297387, 4972.2187531479649, 10142.326190192493, 1012.621416177032, 512.63940465703138, 2.3143394563933557e-05, 2.7811199599886778e-05, 4.5270475247589694e-05, 9.3721328691165031e-05, 4.4757048960255865e-05, 0.00034067378168311287, 9.3019688674096033e-05, 0.00029748422713635892]
pts = [20, 30, 40]
ns = [20, 20, 20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
