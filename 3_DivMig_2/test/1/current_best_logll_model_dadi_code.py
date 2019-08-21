#current best params = [9989.8597311885387, 13254.579666842581, 5179.6582631480323, 16114.516866452037, 4960.964704011526, 10093.046957215538, 1114.7911879070746, 482.75216155154862, 0.00012013273452295402, 1.9328936602124865e-05, 3.5271213084713667e-05, 7.0370582534391241e-05, 4.8628468893893049e-05, 0.00028604144568468854, 9.8687003706937357e-05, 0.00027269742182715564]
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

popt = [9989.8597311885387, 13254.579666842581, 5179.6582631480323, 16114.516866452037, 4960.964704011526, 10093.046957215538, 1114.7911879070746, 482.75216155154862, 0.00012013273452295402, 1.9328936602124865e-05, 3.5271213084713667e-05, 7.0370582534391241e-05, 4.8628468893893049e-05, 0.00028604144568468854, 9.8687003706937357e-05, 0.00027269742182715564]
pts = [20, 30, 40]
ns = [20, 20, 20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
