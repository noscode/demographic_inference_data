import moments
from demographic_model import model_func as func

theta = 20000
popt = [1.0, 2.0, 1.5, 1.0, 0.5, 0.05, 0.1, 0.15, 0.05]

data = func(popt, [20, 20, 20, 20, 20]) * theta
model = data

model.to_file('simulated_data.fs')
print("Data saved")
ll_model = moments.Inference.ll_multinom(model, model)
print('Maximum log composite likelihood: {0}'.format(ll_model))

