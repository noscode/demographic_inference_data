import moments
from demographic_model import model_func as func
import sys
import os
import numpy as np
import multiprocessing as mp
import signal
import random
import Optimize_Functions

id_num = int(sys.argv[1])

def run_one_job(prefix):	
	np.random.seed()
	path_to_dir = os.path.join('results', 'dadi_pipeline', prefix)
	out_file = os.path.join(path_to_dir, 'out.log')
	data = moments.Spectrum.from_file(os.path.join('simulated_data', str(id_num).zfill(2) + '.fs'))
	
	if os.path.exists(path_to_dir):
		if os.path.exists(os.path.join(path_to_dir, 'popt.npy')):
			popt = np.load(os.path.join(path_to_dir, 'popt.npy'))
			model = func(popt, [20])
			ll = moments.Inference.ll_multinom(model, data)
			return int(prefix), ll, popt
	else:
		os.mkdir(path_to_dir)
	
	sys.stdout = open(path_to_dir + '/log_file.log', 'w')
	params = [0.01, 1.0, 0.005, 0.05]
	lower_bound = [0.001, 0.001, 0, 0]
	upper_bound = [100, 100, 5, 5]

	Optimize_Functions.Optimize_Routine(data, [40,50,50], path_to_dir + '/prefix', 'our_model', func, 5, 4, False)
	with open(os.path.join(path_to_dir, 'prefix.our_model.optimized.txt')) as f:
		best_ll = None
		for line in f:
			try:
				ll = float(line.split()[2])
			except:
				continue
			if best_ll is None or best_ll < ll:
				best_ll, best_params = ll, [float(x) for x in line.split()[-1].split(',')]
	popt = np.array(best_params)
	np.save(os.path.join(path_to_dir, 'popt.npy'), popt)
	model = func(popt, [20])
        ll = moments.Inference.ll_multinom(model, data)
	return int(prefix), ll, popt

def worker_init():
        """Graceful way to interrupt all processes by Ctrl+C."""
        # ignore the SIGINT in sub process
        def sig_int(signal_num, frame):
                pass
        signal.signal(signal.SIGINT, sig_int)


pool = mp.Pool(processes=10, initializer=worker_init)
try:
        result = []
#2489,100000000
        map_result = pool.map_async(run_one_job, [str(i) for i in xrange(1, 51)], callback=result.extend)
        while True:
                try:
                        map_result.get(timeout=60)
                        break
                except mp.TimeoutError as ex:
			print map_result
                        pass
                except Exception, e:
                        pool.terminate()
                        sys.stderr.write('ERROR:' + str(e) + '\n')
			os._exit(1)
#	for x in result:
#		print x
	result = sorted(result, key=lambda x: x[1], reverse=True)
	for x in result[:4]:
		print x
	import numpy as np
	a = np.array([x[1] for x in result])
	print np.mean(a), np.std(a)
	print np.mean(a[:50]), np.std(a[:50])
        pool.close()
        pool.join()
except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)


