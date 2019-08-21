import os, sys
import numpy as np

dp_times = []
dp_iters = []
g_times = []
g_iters = []

for i in xrange(1, 51):
	if not os.path.isfile('results/dadi_pipeline/' + str(i) + '/popt.npy'):
		continue
	with open('results/dadi_pipeline/' + str(i) + '/log_file.log') as f:
		for line in f:
			if line.startswith('our_model Analysis Time for Model:'):
				time = line.split()[-2].split(':')
		print i, time
		time = int(time[0]) * 60 * 60 + int(time[1]) * 60 + float(time[2])
	with open('results/dadi_pipeline/' + str(i) + '/prefix.our_model.log.txt') as f:
                last_iter = 0
		for line in f:
			if line.strip() == '' or not line.split()[0].isdigit():
				continue
			last_iter += 1
	dp_times.append(time)
	dp_iters.append(last_iter)
	
	with open('results/gadma/out_0_custom_triangular/' + str(i) + '/GADMA_GA.log') as f:
                for line in f:
			if line.startswith('Iteration '):
				last_iter = int(line.split()[1][1:-1])
			if line.startswith('Total time:'):
				time = line.split()[-2].split(':')
                time = int(time[0]) * 60 * 60 + int(time[1]) * 60 + float(time[2])
	with open('results/gadma/out_0_custom/' + str(i) + '/GADMA_GA.log') as f:
                for line in f:
                        if line.startswith('Iteration '):
                                last_iter = int(line.split()[1][1:-1])
	for fn in os.listdir('results/gadma/out_0_custom/' + str(i) ):
		if fn.startswith('opt'):
			with open('results/gadma/out_0_custom/' + str(i) + '/' + fn) as f:
				cur_it = 0
				for line in f:
					cur_it += 1
	last_iter = last_iter * 10 + cur_it
	g_times.append(time)
	g_iters.append(last_iter)

print dp_times[1], np.mean(dp_times)
print g_times[1], np.mean(g_times)
print dp_iters[1], np.mean(dp_iters)
print g_iters[1], np.mean(g_iters)
