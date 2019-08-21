#current best params = [9845.1209359532404, 93851.252613459263, 4184.8284910254724, 5900.1726280049143, 4952.0078513346307, 11290.438575056833, 1221.4306183828669, 366.31524901795143, 0.0, 0.00028457890089179846, 4.5375404451365484e-05, 0.00011441288541497902, 2.8185550585132386e-05, 0.0, 7.5210908736372508e-05, 2.0696485350407467e-05]
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

popt = [9845.1209359532404, 93851.252613459263, 4184.8284910254724, 5900.1726280049143, 4952.0078513346307, 11290.438575056833, 1221.4306183828669, 366.31524901795143, 0.0, 0.00028457890089179846, 4.5375404451365484e-05, 0.00011441288541497902, 2.8185550585132386e-05, 0.0, 7.5210908736372508e-05, 2.0696485350407467e-05]
pts = [20, 30, 40]
ns = [20, 20, 20]
func_ex = dadi.Numerics.make_extrap_log_func(generated_model)
model =  func_ex(popt, ns, pts)
ll_model = dadi.Inference.ll(model, data)
print('Model log likelihood (LL(model, data)): {0}'.format(ll_model))
