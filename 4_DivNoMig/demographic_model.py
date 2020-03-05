import moments
import numpy as np


def model_func(params, ns):
	nu1, nu234, nu2, nu34, nu3, nu4, T1, T2, T3 = params
	sts = moments.LinearSystem_1D.steady_state_1D(sum(ns))
	fs = moments.Spectrum(sts)

	fs = moments.Manips.split_1D_to_2D(fs, ns[0], sum(ns[1:]))

	fs.integrate(Npop=[nu1, nu234], tf=T1)
	
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], sum(ns[2:]))
	
	fs.integrate(Npop=[nu1, nu2, nu34], tf=T2)

	fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], sum(ns[3:]))
	
	fs.integrate(Npop=[nu1, nu2, nu3, nu4], tf=T3)

	return fs
	
	
