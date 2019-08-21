#current best params = [10000.354374509727, 10016.57555862245, 994.58982019543839, 493.79336255845499, 0.00049796195226443972, 0.00024577194362113066]
import dadi
import numpy as np

def generated_model(params, ns, pts):
	Ns = params[:3]
	Ts = params[3:4]
	Ms = params[4:]
	theta1 = 2.0
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx, theta0=theta1, nu=Ns[0])
	phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
	T = Ts[0]
	after = Ns[1:3]
	growth_func_1 = after[0]
	growth_func_2 = after[1]
	phi = dadi.Integration.two_pops(phi, xx,  T=T,nu1=growth_func_1, nu2=growth_func_2, m12=Ms[0], m21=Ms[1], theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*2)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/2_DivMigr/simulated_data/00.fs')

popt = [10000.354374509727, 10016.57555862245, 994.58982019543839, 493.79336255845499, 0.00049796195226443972, 0.00024577194362113066]
pts = [20, 30, 40]
ns = [20, 20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
