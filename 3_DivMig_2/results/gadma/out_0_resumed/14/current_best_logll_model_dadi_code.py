#current best params = [3121.3506828137083, 118.38763233475149, 2958.1586039585686, 0.021729657495044304, 13611.963344370251, 1935.8475687209079, 7135.4281660068946, 12455.067551026405, 6808.3795050577228, 0.0, 0.0, 0.00020079560134750286, 3.7490515077963783e-05, 0.0, 0.0029267528410053629, 0.00019315565445954418, 0.0]
import dadi
import numpy as np

def generated_model(params, ns, pts):
	Ns = params[:7]
	Ts = params[7:9]
	Ms = params[9:]
	theta1 = 2.0
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx, theta0=theta1, nu=Ns[0])
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	T = Ts[0]
	after = Ns[1:3]
	growth_func_1 = after[0]
	growth_func_2 = after[1]
	phi = dadi.Integration.two_pops(phi, xx,  T=T,nu1=growth_func_1, nu2=growth_func_2, m12=Ms[0], m21=Ms[1], theta0=theta1)
	before = after
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	before.append((1 - Ns[3]) * before[-1])
	before[-2] *= Ns[3]
	T = Ts[1]
	after = Ns[4:7]
	growth_func_1 = lambda t: before[0] + (after[0] - before[0]) * (t / T)
	growth_func_2 = after[1]
	growth_func_3 = lambda t: before[2] * (after[2] / before[2]) ** (t / T)
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=Ms[2], m13=Ms[3], m21=Ms[4], m23=Ms[5], m31=Ms[6], m32=Ms[7], theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [3121.3506828137083, 118.38763233475149, 2958.1586039585686, 0.021729657495044304, 13611.963344370251, 1935.8475687209079, 7135.4281660068946, 12455.067551026405, 6808.3795050577228, 0.0, 0.0, 0.00020079560134750286, 3.7490515077963783e-05, 0.0, 0.0029267528410053629, 0.00019315565445954418, 0.0]
pts = [20, 30, 40]
ns = [20, 20, 20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
