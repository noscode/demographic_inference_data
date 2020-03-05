import moments
import numpy as np


def model_func(params, ns):
	nu1, nu2, nu3, nu4, nu5, T1, T2, T3, T4 = params
	sts = moments.LinearSystem_1D.steady_state_1D(sum(ns))
	fs = moments.Spectrum(sts)

	fs = moments.Manips.split_1D_to_2D(fs, ns[0], sum(ns[1:]))

	fs.integrate(Npop=[nu1, nu2], tf=T1)
	
	fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], sum(ns[2:]))
	
	fs.integrate(Npop=[nu1, nu2, nu3], tf=T2)

	fs = moments.Manips.split_3D_to_4D_3(fs, ns[2], sum(ns[3:]))
	
	fs.integrate(Npop=[nu1, nu2, nu3, nu4], tf=T3)

	fs = moments.Manips.split_4D_to_5D_4(fs, ns[3], ns[4])

        fs.integrate(Npop=[nu1, nu2, nu3, nu4, nu5], tf=T4)

	return fs
	
	
