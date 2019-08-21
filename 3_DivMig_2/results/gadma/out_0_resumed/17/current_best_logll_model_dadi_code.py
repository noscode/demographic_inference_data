#current best params = [10126.092521942077, 0.79617241371674641, 7365.7859840467172, 3648.4156753471652, 0.57812963138303286, 15956.680256129115, 4863.4845393689393, 12431.054161907363, 5564.9262054603114, 1217.3504614242142, 0.00067645276931379032, 0.00098754776122489006, 2.9595426134274518e-05, 6.8460559673336394e-05, 3.552270116437394e-05, 0.00076506501994233532, 8.7532955250038197e-05, 0.00079792922809833759]
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
	growth_func_2 = lambda t: before[1] + (after[1] - before[1]) * (t / T)
	phi = dadi.Integration.two_pops(phi, xx,  T=T,nu1=growth_func_1, nu2=growth_func_2, m12=Ms[0], m21=Ms[1], theta0=theta1)
	before = after
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	before.append((1 - Ns[4]) * before[-1])
	before[-2] *= Ns[4]
	T = Ts[1]
	after = Ns[5:8]
	growth_func_1 = after[0]
	growth_func_2 = lambda t: before[1] + (after[1] - before[1]) * (t / T)
	growth_func_3 = lambda t: before[2] * (after[2] / before[2]) ** (t / T)
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=Ms[2], m13=Ms[3], m21=Ms[4], m23=Ms[5], m31=Ms[6], m32=Ms[7], theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [10126.092521942077, 0.79617241371674641, 7365.7859840467172, 3648.4156753471652, 0.57812963138303286, 15956.680256129115, 4863.4845393689393, 12431.054161907363, 5564.9262054603114, 1217.3504614242142, 0.00067645276931379032, 0.00098754776122489006, 2.9595426134274518e-05, 6.8460559673336394e-05, 3.552270116437394e-05, 0.00076506501994233532, 8.7532955250038197e-05, 0.00079792922809833759]
pts = [20, 30, 40]
ns = [20, 20, 20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
