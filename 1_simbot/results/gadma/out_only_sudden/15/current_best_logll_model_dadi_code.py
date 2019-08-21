#current best params = [10092.317871615642, 427.63563427366512, 10476.963411925934, 231.28118719032122, 454.74633388214215]
import dadi
import numpy as np

def generated_model(params, ns, pts):
	Ns = params[:3]
	Ts = params[3:5]
	Ms = params[5:]
	theta1 = 2.0
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx, theta0=theta1, nu=Ns[0])
	T = Ts[0]
	after = Ns[1:2]
	growth_func = after[0]
	phi = dadi.Integration.one_pop(phi, xx, nu=growth_func,T=T, theta0=theta1)
	T = Ts[1]
	after = Ns[2:3]
	growth_func = after[0]
	phi = dadi.Integration.one_pop(phi, xx, nu=growth_func,T=T, theta0=theta1)
	sfs = dadi.Spectrum.from_phi(phi, ns, [xx]*1)
	return sfs
data = dadi.Spectrum.from_file('/mnt/ssd1/enoskova/simulations/1_simbot/SimBot.fs')

popt = [10092.317871615642, 427.63563427366512, 10476.963411925934, 231.28118719032122, 454.74633388214215]
pts = [20, 30, 40]
ns = [20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
