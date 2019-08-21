#current best params = [16437.57576974465, 0.26185605111073817, 201.96913916934525, 10168.031482352897, 0.03693520965488966, 17323.020809559282, 3193.496889172457, 5309.5946340396094, 21269.597348858755, 11102.612139861732, 0.00020357181484132169, 6.8901091493163223e-06, 0.00027385160540629571, 4.0768157965384667e-07, 0.0, 0.00139651274355663, 0.00016313114930937859, 0.00077561089090817442]
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
	growth_func_2 = lambda t: before[1] * (after[1] / before[1]) ** (t / T)
	growth_func_3 = after[2]
	phi = dadi.Integration.three_pops(phi, xx,  T=T, nu1=growth_func_1, nu2=growth_func_2, nu3=growth_func_3, m12=Ms[2], m13=Ms[3], m21=Ms[4], m23=Ms[5], m31=Ms[6], m32=Ms[7], theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*3)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/3_DivMig_2/simulated_data/00.fs')

popt = [16437.57576974465, 0.26185605111073817, 201.96913916934525, 10168.031482352897, 0.03693520965488966, 17323.020809559282, 3193.496889172457, 5309.5946340396094, 21269.597348858755, 11102.612139861732, 0.00020357181484132169, 6.8901091493163223e-06, 0.00027385160540629571, 4.0768157965384667e-07, 0.0, 0.00139651274355663, 0.00016313114930937859, 0.00077561089090817442]
pts = [20, 30, 40]
ns = [20, 20, 20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
