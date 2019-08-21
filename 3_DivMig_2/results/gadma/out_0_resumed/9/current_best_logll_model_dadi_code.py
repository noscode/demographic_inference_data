#current best params = [9997.7715473011449, 14933.714582964159, 4803.6442891832949, 0.99791248490025131, 15028.771811051089, 5163.6246455597984, 10079.640575784822, 1005.5757848057365, 508.09677841654081, 0.0, 4.4016925212442216e-05, 4.7887023825785343e-05, 0.00010086760240868096, 4.2704766759326242e-05, 0.00031811426407407818, 8.9486321161522117e-05, 0.00031012211652042154]
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
	growth_func_2 = lambda t: before[1] * (after[1] / before[1]) ** (t / T)
	growth_func_3 = after[2]
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=Ms[2], m13=Ms[3], m21=Ms[4], m23=Ms[5], m31=Ms[6], m32=Ms[7], theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [9997.7715473011449, 14933.714582964159, 4803.6442891832949, 0.99791248490025131, 15028.771811051089, 5163.6246455597984, 10079.640575784822, 1005.5757848057365, 508.09677841654081, 0.0, 4.4016925212442216e-05, 4.7887023825785343e-05, 0.00010086760240868096, 4.2704766759326242e-05, 0.00031811426407407818, 8.9486321161522117e-05, 0.00031012211652042154]
pts = [20, 30, 40]
ns = [20, 20, 20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
