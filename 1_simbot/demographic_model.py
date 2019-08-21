import moments

import numpy as np

def model_func((nuB, nuF, tB, tF), ns):
	sts = moments.LinearSystem_1D.steady_state_1D(ns[0])
	fs = moments.Spectrum(sts)
	fs.integrate([nuB], tB)
	fs.integrate([nuF], tF)
	return fs
