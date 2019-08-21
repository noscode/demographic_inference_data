#current best params = [10004.884123469577, 0.12543216486423883, 19813.917129352943, 3334.2999807370602, 0.57972375655833963, 13616.278339394097, 4738.9925674570304, 11643.55163968436, 246.56625229672946, 1371.3923429781776, 0.00018903734878048662, 0.00020029459159407549, 4.8701771299017758e-05, 6.0888272830147306e-05, 3.8967493862144323e-05, 0.00079031037844855723, 9.727953123604356e-05, 0.00080210451906456942]
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
	growth_func_1 = lambda t: before[0] + (after[0] - before[0]) * (t / T)
	growth_func_2 = after[1]
	phi = dadi.Integration.two_pops(phi, xx,  T=T,nu1=growth_func_1, nu2=growth_func_2, m12=Ms[0], m21=Ms[1], theta0=theta1)
	before = after
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)
	before.append((1 - Ns[4]) * before[-1])
	before[-2] *= Ns[4]
	T = Ts[1]
	after = Ns[5:8]
	growth_func_1 = lambda t: before[0] * (after[0] / before[0]) ** (t / T)
	growth_func_2 = lambda t: before[1] + (after[1] - before[1]) * (t / T)
	growth_func_3 = lambda t: before[2] * (after[2] / before[2]) ** (t / T)
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=Ms[2], m13=Ms[3], m21=Ms[4], m23=Ms[5], m31=Ms[6], m32=Ms[7], theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [10004.884123469577, 0.12543216486423883, 19813.917129352943, 3334.2999807370602, 0.57972375655833963, 13616.278339394097, 4738.9925674570304, 11643.55163968436, 246.56625229672946, 1371.3923429781776, 0.00018903734878048662, 0.00020029459159407549, 4.8701771299017758e-05, 6.0888272830147306e-05, 3.8967493862144323e-05, 0.00079031037844855723, 9.727953123604356e-05, 0.00080210451906456942]
pts = [20, 30, 40]
ns = [20, 20, 20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
