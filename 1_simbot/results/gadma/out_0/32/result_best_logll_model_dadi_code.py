#current best params = [10190.31156602354, 269.27375281575729, 28598.086794846738, 3.0876470667529298e-10, 756.47739910339965]
import dadi
import numpy as np

def generated_model(params, ns, pts):
	Ns = params[:3]
	Ts = params[3:5]
	Ms = params[5:]
	theta1 = 2.0
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx, theta0=theta1, nu=Ns[0])
	before = [Ns[0]]
	T = Ts[0]
	after = Ns[1:2]
	growth_func = lambda t: before[0] * (after[0] / before[0]) ** (t / T)
	phi = dadi.Integration.one_pop(phi, xx, nu=growth_func,T=T, theta0=theta1)
	before = after
	T = Ts[1]
	after = Ns[2:3]
	growth_func = lambda t: before[0] * (after[0] / before[0]) ** (t / T)
	phi = dadi.Integration.one_pop(phi, xx, nu=growth_func,T=T, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*1)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/1_simbot/simulated_data/00.fs')

popt = [10190.31156602354, 269.27375281575729, 28598.086794846738, 3.0876470667529298e-10, 756.47739910339965]
pts = [20, 30, 40]
ns = [20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))