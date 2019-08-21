#current best params = [10060.2253535073, 0.78119533145150966, 34778.113098094385, 10546.024389682447, 0.50945947763997879, 3348.5808421657448, 4285.5350216377146, 12163.809500597694, 1112.3323248553168, 341.99366805514808, 5.1461772907270248e-06, 0.0, 6.7054119021242473e-05, 0.00014151304977544829, 6.7947610583648195e-05, 0.00011912166643842915, 0.00013839309354127111, 0.00017893863068120855]
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
	growth_func_1 = lambda t: before[0] * (after[0] / before[0]) ** (t / T)
	growth_func_2 = lambda t: before[1] + (after[1] - before[1]) * (t / T)
	phi = dadi.Integration.two_pops(phi, xx,  T=T,nu1=growth_func_1, nu2=growth_func_2, m12=Ms[0], m21=Ms[1], theta0=theta1)
	before = after
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	before.append((1 - Ns[4]) * before[-1])
	before[-2] *= Ns[4]
	T = Ts[1]
	after = Ns[5:8]
	growth_func_1 = lambda t: before[0] + (after[0] - before[0]) * (t / T)
	growth_func_2 = after[1]
	growth_func_3 = lambda t: before[2] + (after[2] - before[2]) * (t / T)
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=Ms[2], m13=Ms[3], m21=Ms[4], m23=Ms[5], m31=Ms[6], m32=Ms[7], theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [10060.2253535073, 0.78119533145150966, 34778.113098094385, 10546.024389682447, 0.50945947763997879, 3348.5808421657448, 4285.5350216377146, 12163.809500597694, 1112.3323248553168, 341.99366805514808, 5.1461772907270248e-06, 0.0, 6.7054119021242473e-05, 0.00014151304977544829, 6.7947610583648195e-05, 0.00011912166643842915, 0.00013839309354127111, 0.00017893863068120855]
pts = [20, 30, 40]
ns = [20, 20, 20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
