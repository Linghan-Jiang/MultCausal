import os
import random
from multiprocessing import Pool

def gsmr(bfile='data/arabi', out='results/gsmr/test', gsmr_d='data/gsmr', gwas_thresh=5e-08, n=20, test=None):
	if not os.path.isdir(out):
		os.makedirs(out)
	exposure_d = '{}/exposure'.format(gsmr_d)
	outcome_d = '{}/outcome'.format(gsmr_d)
	exposure_f = '{}/gsmr_exposure.txt'.format(gsmr_d)
	outcome_f = '{}/gsmr_outcome.txt'.format(gsmr_d)
	tmp_d = '{}/tmp_{}'.format(gsmr_d, str(random.random())[-10:])
	print(tmp_d)
	os.makedirs(tmp_d)
	order = 'split -l 350 -a 3 {} {}/'.format(exposure_f, tmp_d)
	os.system(order)
	lst = os.listdir(tmp_d)
	if test:
		n = 1
		lst = lst[:1]
	pool = Pool(n) 
	for i in lst:
		out_f = '{}_{}'.format(out, i)
		f1 = '{}/{}'.format(tmp_d, i)
		order = 'gcta64 --bfile {} --gsmr-file {} {} --gsmr-direction 0 --gwas-thresh {} --out {}'.format(bfile, f1, outcome_f, gwas_thresh, out_f) 
		#os.system(order)
		pool.apply_async(os.system, (order,))
	pool.close()
	pool.join()
	os.system('rm -rf {}'.format(tmp_d))
if __name__ == '__main__':
	for i in range(7, 1, -1):
		out = 'data/gsmr/gsmr_5e-0{}/test'.format(i)
		gwas_thresh = '5e-0{}'.format(i)	
		gsmr(out=out, gwas_thresh=gwas_thresh)
