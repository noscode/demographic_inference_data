import moments
from demographic_model import model_func as func
import sys
import os
import numpy as np
import multiprocessing as mp
import signal
import random

id_num = int(sys.argv[1])

def run_one_job(prefix):	
	np.random.seed()
	path_to_dir = os.path.join('results', 'moments', str(id_num) + '_random', prefix)
	out_file = os.path.join(path_to_dir, 'out.log')
	data = moments.Spectrum.from_file(os.path.join('simulated_data', str(id_num).zfill(2) + '.fs'))
	
	if os.path.exists(path_to_dir):
		if os.path.exists(os.path.join(path_to_dir, 'popt.npy')):
			popt = np.load(os.path.join(path_to_dir, 'popt.npy'))
			model = func(popt, [20, 20, 20])
			ll = moments.Inference.ll_multinom(model, data)
			return int(prefix), ll, popt
	else:
		os.mkdir(path_to_dir)
	

	params = [1.5, 0.5, 1.0, 0.5, 1.0, 3.0, 0.1, 0.05]
	lower_bound = [0.001, 0.001, 1e-2, 0,0,0, 0, 0]
	upper_bound = [100, 100, 100, 10, 10, 10, 5, 5]

	p0 = params	
	p0 = moments.Misc.perturb_params(params, fold=3, upper_bound=upper_bound, lower_bound=lower_bound)
	p0 = [np.exp(np.random.triangular(np.log(lower_bound[i] + 1e-100), np.log(1.0), np.log(upper_bound[i]))) for i in xrange(len(upper_bound))]	
	#print('Beginning optimization ************************************************')
	popt = moments.Inference.optimize_powell(p0, data, func,
	                                   lower_bound=lower_bound,
	                                   upper_bound=upper_bound,
	                                   verbose=1, maxiter=100, output_file=out_file)
	# The verbose argument controls how often progress of the optimizer should be
	# printed. It's useful to keep track of optimization process.
	#print('Finished optimization **************************************************')
	popt = np.array(popt)
	np.save(os.path.join(path_to_dir, 'popt.npy'), popt)
	model = func(popt, [20, 20, 20])
        ll = moments.Inference.ll_multinom(model, data)
	return int(prefix), ll, popt

def worker_init():
        """Graceful way to interrupt all processes by Ctrl+C."""
        # ignore the SIGINT in sub process
        def sig_int(signal_num, frame):
                pass
        signal.signal(signal.SIGINT, sig_int)

if not os.path.exists(os.path.join('results', 'moments', str(id_num) + '_random')):
	os.mkdir(os.path.join(os.path.join('results', 'moments', str(id_num) + '_random')))

pool = mp.Pool(processes=10, initializer=worker_init)
try:
        result = []
        map_result = pool.map_async(run_one_job, [str(i) for i in xrange(1,251)], callback=result.extend)
        while True:
                try:
                        map_result.get(timeout=60)
                        break
                except mp.TimeoutError as ex:
			print map_result
                        pass
                except Exception, e:
                        pool.terminate()
			print e
                        sys.stderr.write('ERROR:' + str(e) + '\n')
			raise e
			os._exit(1)
#	for x in result:
#		print x
	result = sorted(result, key=lambda x: x[1], reverse=True)
	for x in result[:4]:
		print x
	lls = []
	for x in result:
		lls.append(x[1])
	import numpy as np
	lls = np.array(lls)
	print np.mean(lls), np.std(lls)
        pool.close()
        pool.join()
except KeyboardInterrupt:
        pool.terminate()
        sys.exit(1)


