import moments
import numpy as np
import time
from demographic_model import model_func

# Parameters from Jouganous et al. (2017) paper (Table 4):
N_A = 11293
N_Af = 24486
N_B = 3034
N_Eu0 = 2587
r_Eu = 0.17e-2  # it is percent in table
N_As0 = 958
r_As = 0.30e-2  # it is percent in table
m_Af_B = 15.6e-5
m_Af_Eu = 1.00e-5
m_Af_As = 0.48e-5
m_Eu_As = 3.99e-5
T_Af = 349e3  # kya
T_B = 121e3  # kya
T_Eu_As = 44e3  # kya

# Time for one generation and mutation rate in paper
t_gen = 29  # years
mu = 1.44e-8

# translate these values to the one from model
# 1. Translate time to time of intervals and
T_Af -= T_B
T_B -= T_Eu_As
# Translate it to generations
T_g_Af = T_Af / t_gen
T_g_B = T_B / t_gen
T_g_Eu_As = T_Eu_As / t_gen
# 2. Get final sizes from r
# (1 + r) ** t = N_final / N_init
# N_final = N_init * ((1 + r) ** t)
# t is in generations!
N_Eu = N_Eu0 * ((1 + r_Eu) ** T_g_Eu_As)
N_As = N_As0 * ((1 + r_As) ** T_g_Eu_As)
print(N_Eu, N_As)
# 3. Get relative sizes
Npop = np.array([N_Af, N_B, N_Eu0, N_Eu, N_As0, N_As])
nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs = Npop / N_A
# 4. Translate migrations by * 2 * N_A
m = np.array([m_Af_B, m_Af_Eu, m_Af_As, m_Eu_As])
mAfB, mAfEu, mAfAs, mEuAs = m * 2 * N_A
# 5. Translate time from generations by / (2 * N_A)h
TAf = T_g_Af / (2 * N_A)
TB = T_g_B / (2 * N_A)
TEuAs = T_g_Eu_As / (2 * N_A)
# Form params
params = (nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs,
          mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs)

print(params)

# Load data
data = moments.Spectrum.from_file("fs_data_40.fs")

ns = (40, 40, 40)#data.sample_sizes

# Draw model
#from matplotlib import pyplot as plt
#model = moments.ModelPlot.generate_model(model_func, params, (10, 10, 10, 10))
#moments.ModelPlot.plot_model(model, 
#	save_file='model.png',
#	fig_title='',
#	pop_labels=['YRI', 'CEU', 'CHB', "JPT"],
#	nref=N_A,
#	draw_scale=True,
#	gen_time=0.029,
#	gen_time_units="Thousand years",
#	reverse_timeline=True)
#plt.show()

# Simulate
t1 = time.time()
model = model_func(params, ns)
t2 = time.time()
print(f"Time of simulation is {(t2 - t1):2f} seconds.")

# Calculate ll
ll = moments.Inference.ll_multinom(model, data)
print(f"Log likelihood of model is equal to {ll:2f}.")

#lower_bound = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 0, 0, 0, 0, 0, 0, 0, 0, 0]
#upper_bound = [100, 100, 100, 100, 100, 100, 100, 100, 10, 10, 10, 10, 10, 5, 5, 5, 5]
#
#assert len(params) == len(lower_bound)
#assert len(params) == len(upper_bound)
#
#popt = moments.Inference.optimize_log(model_func=model_func, data=data, p0=params, lower_bound=lower_bound, upper_bound=upper_bound, verbose=1)
