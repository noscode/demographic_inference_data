import os, sys
import numpy as np

dirname = sys.argv[1]

d = {}
start_times = {}
with open(os.path.join(dirname, 'GADMA.log')) as f:
	for line in f:
		if line.startswith('['):
			time = line.strip()[1:-1].split(':')
			time = int(time[0]) * 60.0 * 60 + int(time[1]) * 60 + int(time[2])
		if line.strip() != '' and line.split()[0] == 'Model' and  line.split()[1].isdigit():
			num = int(line.split()[1])
			info = line.strip().split()[-1]
			if num not in start_times:
				start_times[num] = time
			if info == 'f' and num not in d:
				d[num] = time - start_times[num]
iters =[]
for i in xrange(1, 11):
        with open(os.path.join(dirname,  str(i), 'GADMA_GA.log')) as f:
                for line in f:
                        if line.startswith('Iteration '):
                                last_iter = int(line.split()[1][1:-1])
	cur_it = 0
        for fn in os.listdir(os.path.join(dirname, str(i))):
                if fn.startswith('opt'):
                        with open(os.path.join(dirname, str(i),  fn)) as f:
                                cur_it = 0
                                for line in f:
                                        cur_it += 1
        last_iter = last_iter * 10 + cur_it
	iters.append(last_iter)
import numpy as np
a = np.array([d[x] for x in d])
print 'time:', np.mean(a)
print 'iter:', np.mean(iters)
