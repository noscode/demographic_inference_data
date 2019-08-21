import numpy as np
import sys
import os
import multiprocessing as mp
import signal
import random

moments_model_filename = 'demographic_model.py'
par_ident = {1: 'n,n,t,t', 2: 'n,n,n,m,m,t,t'}
struct = {1: '3', 2: '1,1'}

def create_gadma_params(id_num, custom, sudden_only):
	suffix = ''
	if custom:
		suffix += '_custom'
	if sudden_only:
		suffix += '_sudden_only'
	with open(os.path.join('results','gadma', 'params' + str(id_num) + suffix), 'w') as f:
		f.write('Output directory: results/gadma/out_' + str(id_num) + suffix + '\n')
		f.write('Input file: simulated_data/' + str(id_num).zfill(2) + '.fs' + '\n')
		f.write('Use moments or dadi: moments' + '\n')
		f.write('Theta0: 2.0\n')
		if custom:
			f.write('Custom filename: ' + moments_model_filename + '\n')
			f.write('Parameter identifiers: ' + par_ident[1] + '\n')
		else:
			f.write('Initial structure: ' + struct[1] + '\n')
		if sudden_only:
			f.write('Only sudden: True' + '\n')
		f.write('Number of processes: 10' + '\n')
		f.write('Number of repeats: 50' + '\n')
		f.write('Silence: True\n')
		f.write('min_N: 1e-3')
	return os.path.join('results','gadma', 'params' + str(id_num) + suffix + '\n')
cmds = []
for i in xrange(1):
	params_file = create_gadma_params(i, False, False)
	cmd_gadma = 'gadma -p ' + params_file
	cmds.append(cmd_gadma)
	
	params_file = create_gadma_params(i, True, False)
	cmd_gadma = 'gadma -p ' + params_file
	cmds.append(cmd_gadma)

	params_file = create_gadma_params(i, False, True)
	cmd_gadma = 'gadma -p ' + params_file
	cmds.append(cmd_gadma)
	
	cmd_moments = 'python run_moments.py ' + str(i)
	cmds.append(cmd_moments)

def run_one_job(cmd):
	np.random.seed()
	print cmd
	import subprocess
	process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
	output, error = process.communicate()


def worker_init():
        """Graceful way to interrupt all processes by Ctrl+C."""
        # ignore the SIGINT in sub process
        def sig_int(signal_num, frame):
                pass
        signal.signal(signal.SIGINT, sig_int)

pool = mp.Pool(processes=2, initializer=worker_init)
try:
	result = []
	map_result = pool.map_async(run_one_job, cmds, callback=result.extend)
	while True:
		try:
			map_result.get(timeout=1)
			break
		except mp.TimeoutError as ex:
			pass
		except Exception, e:
			pool.terminate()
			sys.stderr.write(str(e) + '\n')
			os._exit(1)
	pool.close()
	pool.join()
except KeyboardInterrupt:
	pool.terminate()
	sys.exit(1)

