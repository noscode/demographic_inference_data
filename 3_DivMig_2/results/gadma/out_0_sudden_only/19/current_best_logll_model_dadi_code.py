#current best params = [9966.5266487159843, 12978.326261081078, 5348.5381388481173, 15980.712278915997, 4957.6719563598626, 9921.9457897074517, 1178.7771760806104, 473.82362424774078, 0.00016027807123787896, 1.4782847613878896e-05, 3.1027064033337124e-05, 6.2509173140633972e-05, 5.0420882821972815e-05, 0.00027415803328018523, 0.00010332733015642413, 0.00026989185296157233]
import dadi
import numpy as np

def generated_model(params, ns, pts):
	Ns = params[:6]
	Ts = params[6:8]
	Ms = params[8:]
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
	T = Ts[1]
	after = Ns[3:6]
	growth_func_1 = after[0]
	growth_func_2 = after[1]
	growth_func_3 = after[2]
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=Ms[2], m13=Ms[3], m21=Ms[4], m23=Ms[5], m31=Ms[6], m32=Ms[7], theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [9966.5266487159843, 12978.326261081078, 5348.5381388481173, 15980.712278915997, 4957.6719563598626, 9921.9457897074517, 1178.7771760806104, 473.82362424774078, 0.00016027807123787896, 1.4782847613878896e-05, 3.1027064033337124e-05, 6.2509173140633972e-05, 5.0420882821972815e-05, 0.00027415803328018523, 0.00010332733015642413, 0.00026989185296157233]
pts = [20, 30, 40]
ns = [20, 20, 20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
