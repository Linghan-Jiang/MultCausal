import os
import scipy
import random
import subprocess
import numpy as np
import pandas as pd
from utils import *
from multiprocessing import Pool
def worker(ind, f1, f2, out_f, step=1000): # f save in 'scipy.sparese.npz' format
		m1 = scipy.sparse.load_npz(f1) 
		m2 = scipy.sparse.load_npz(f2)
		m2 = csr_matrix(m2[ind*step:(ind+1)*step, :])
		if m1.shape[1] != m2.shape[0]:
			raise('dimention incompatible!')
		m = m1.dot(m2)	
		scipy.sparse.save_npz(out_f, m)
def dd_dot(d1, d2, out, parallel=True):
	print('Starting dd_dot ...')
	print('d1: ', d1)
	print('d2: ', d2)
	d1_lst = sorted(os.listdir(d1)) # sg_d
	d2_lst = sorted(os.listdir(d2)) # si_d
	if (d1_lst == d2_lst):
		print('d1 and d2 have the same order.')
	else:
		print('len(d1):', len(d1))
		print('len(d2):', len(d2))
		raise("ERROR: d1 doesn't match with d2!")
	print('ig d:', out)
	out_lst = os.listdir(out)
	out_lst = [i.split('.')[0] for i in out_lst]
	if len(out_lst) != len(d1_lst):
		pool = Pool(10)
		for ind, (d1_f, d2_f) in enumerate(zip(d1_lst, d2_lst)):
			if d1_f not in out_lst:
				common_f = d1_f
				out_f = '{}/{}'.format(out, common_f)
				d1_f = '{}/{}'.format(d1, common_f)
				d2_f = '{}/{}'.format(d2, common_f)
	#			print('into loop')
				if parallel:
					pool.apply_async(worker,(ind, d1_f, d2_f, out_f,))
				else:
					worker(ind, d1_f, d2_f, out_f)
		pool.close()
		pool.join()
		go_on = True
	else:
		go_on = False
		print('dd mult done, results save to {}.'.format(out))
	return go_on 
def mult_casual(pre, out_pre, gg_f=None, sg_pre=None, reg=None):
	gdic_f = '{}_gdic.txt'.format(pre)
	pheno_f = '{}.pheno'.format(pre)
	if sg_pre:
		sg_d = '{}_sg'.format(sg_pre)
	else:
		sg_d = '{}_sg'.format(out_pre)
	is_d = '{}_is'.format(out_pre)
	ig_d = '{}_ig'.format(out_pre)
	out_d, tmp = os.path.split(out_pre)
	if not os.path.isdir(out_d):
		os.makedirs(out_d)
	geno_f = '{}.is.geno'.format(out_pre)
	# main
	if not reg:
		go_on = dd_dot(is_d, sg_d, ig_d, parallel=False)
		if go_on:
			dd_dot(is_d, sg_d, ig_d)
			print('complting igmat ...')
		ig = os.listdir(ig_d)
		igmat = scipy.sparse.load_npz(os.path.join(ig_d, ig[0]))
		for i in ig[1:]:
			igmat += scipy.sparse.load_npz(os.path.join(ig_d, i))
		if gg_f:
			ggmat = scipy.sparse.load_npz(gg_f)
			igmat = igmat.dot(ggmat)
			_, gg_label = os.path.split(gg_f)
			gg_label = '.'.join(gg_label.split('.')[:-2])
			geno_f = '{}.isg.{}.geno'.format(out_pre, gg_label)
			tmp = '{}.isg.{}'.format(tmp, gg_label)
		igmat = igmat.toarray()
		igmat = rescale(igmat)
		gimat = igmat.transpose()
		gimat = pd.DataFrame(gimat)
		print('the final shape is :', gimat.shape)
		lmat = left(gdic_f)
		omat = pd.concat([lmat, gimat], ignore_index=True, axis=1)
		print('outputing geno to {}.'.format(geno_f))
		omat.to_csv(geno_f, header=None, index=False)
	order = 'gemma -g {} -p {} -lm 1 -outdir {} -o {}'.format(geno_f, pheno_f, out_d, tmp); print(order); os.system(order)
def mult_input(pre, out_pre, sg_pre=None):
	si_f = '{}.geno'.format(pre)
	gdic_f = '{}_gdic.txt'.format(pre)
	sdic_f = '{}_sdic.txt'.format(pre)
	is_d = '{}_is'.format(out_pre) 
	si_d = '{}_si'.format(out_pre) 
	ig_d = '{}_ig'.format(out_pre)
	if sg_pre:
		sg_d = '{}_sg'.format(sg_pre)
		sg_f = '{}.sg.txt'.format(sg_pre)
	else:
		sg_d = '{}_sg'.format(out_pre)
		sg_f = '{}.sg.txt'.format(out_pre)
	if not os.path.isdir(is_d):
		os.makedirs(is_d)
	if not os.path.isdir(si_d):
		os.makedirs(si_d)
	if not os.path.isdir(ig_d):
		os.makedirs(ig_d)
	if not os.path.isdir(sg_d):
		os.makedirs(sg_d)
	bimbam_2_is(si_f, is_d=is_d, si_d=si_d)	
	sg(sg_f, si_d=si_d, sg_d=sg_d, gdic_f=gdic_f, sdic_f=sdic_f)
def main(sg_pre, gg_f, out_pre): 
#	r = mult_input('data/mult/arabi', out_pre='data/mult/20kb', sg_pre='data/mult/arabi')
#	r = mult_casual('data/mult/arabi', out_pre='data/mult/20kb', sg_pre='data/mult/arabi', gg_f='data/ppi_mult/atpin.2.no-weight.ppi.txt.npz')
#	r = mult_casual('data/mult/arabi', out_pre='data/mult/20kb', sg_pre='data/mult/arabi', gg_f='data/ppi_mult/atpin.2.weight.ppi.txt.npz')
	_, tmp_pre = os.path.split(sg_pre)	
	_, tmp_pre_2 = os.path.split(gg_f)
	tmp_pre_2 = tmp_pre_2[:-4]
	print(tmp_pre)
	print(tmp_pre_2)
	raise
	r = mult_input('data/mult/arabi', out_pre='data/mult/{}'.format(tmp_pre), sg_pre=sg_pre)
	r = mult_casual('data/mult/arabi', out_pre='data/mult/{}'.format(tmp_pre), sg_pre=sg_pre)
	r = mult_casual('data/mult/arabi', out_pre='data/mult/{}.{}'.format(tmp_pre, tmp_pre_2), sg_pre=sg_pre, gg_f=gg_f)
		 
if __name__ == '__main__':
# 	test
#	mult_input('./test/test')
#	r = mult_casual('./test/test', 'test', reg=True)
#	r = mult_casual('./test/test', 'test', gg_f='test/test.ppi.txt.npz')
	print('Begin.')
	sg_d = 'data/mateqtl'
	ppi_d = 'data/ppi_mult'
	for i in os.listdir(sg_d):
		for j in os.listdir(ppi_d):
			if i.endswith('sg.txt'):
				i = '{}/{}'.format(sg_d, i[:-7])
				j = '{}/{}'.format(ppi_d, j)
				print(i)
				print(j)
				raise
				#main(sg_pre=i, gg_f=j)
