import os
import numpy as np
from numpy.random import *
import pandas as pd
import networkx as nx
import warnings
#from scipy.stats import 
from matplotlib import pyplot as plt
from multiprocessing import Pool
warnings.filterwarnings('ignore')
def f1():
	f = '../data/GSE54680_10C_quant.plink.exp'
	d = pd.read_csv(f, sep=' ', index_col=[0,1])
	causal = np.random.choice(range(1,23447), 10)
	r =  d.iloc[:,causal]
	r.to_csv('')

def causal_snp(n=10):
	f = '../data/arabi.bim'	
	d = pd.read_csv(f,sep='\t',header=None)
	for j in range(n):
		r = np.random.choice(d[1], 10)
		with open('causal_{}.snplist'.format(j), 'w+') as f:
			for i in r:
				f.write(i+'\n')
	
	
	print('into gene_expr')
	bxy = args['bxy']
	bzx = args['bzx']
	rzx2 = args['rzx2']
	rxy2 = args['rxy2']
	c = args['c']
	rc2 = args['rc2']
	if causal_lst:
		m = 1200
#		n = np.random.randint(30, 100)
		n = 600
		r1 = np.random.choice(tmp['SNP'], n, replace=False)
		p = np.array(tmp[tmp['SNP'].isin(r1)]['MAF'])
		w = binomial(1, p)
		r2 = np.random.choice(tmp['SNP'], m-n, replace=False)
		p2 = np.array(tmp[tmp['SNP'].isin(r2)]['MAF'])
		w2 = binomial(1, p2)

		with open('tmp/w_{}'.format(ind), 'w+') as f:
			for i in w:
				f.write(str(i)+'\n')
		with open('tmp/w2_{}'.format(ind), 'w+') as f:
			for i in w2:
				f.write(str(i)+'\n')
		with open('tmp/p_{}'.format(ind), 'w+') as f:
#				f.write(' '.join(str(i))+'\n')
			for i in p:
				f.write(str(i)+'\n')
		with open('tmp/r1_{}'.format(ind), 'w+') as f:
			for i in r1:
				f.write(str(i)+'\n')	
		with open('tmp/r2_{}'.format(ind), 'w+') as f:
			for i in r2:
				f.write(str(i)+'\n')	
	if causal_files:
		r1 = np.array([i.strip() for i in open(causal_files[0]).readlines()])
		r2 = np.array([i.strip() for i in open(causal_files[1]).readlines()])
		p = np.array([float(i.strip()) for i in open(causal_files[2]).readlines()])
		w = np.array([int(i.strip()) for i in open(causal_files[3]).readlines()])
		w2 = np.array([int(i.strip()) for i in open(causal_files[4]).readlines()])
		n = len(r1)
		m = len(r2)
		print(r1, r2, p, w, w2)
#	print(w.shape)
#	print(p.shape)		
	z = (w-2*p)/(np.sqrt(2*p*(1-p)))
	z_bzx = z*bzx
	w_bzx = w*bzx
	var_z_bzx = np.nanvar(z_bzx)
	var_w_bzx = np.nanvar(w_bzx)
	if not rc2:
		rc2 = uniform(0, 0.5)
	sigmac_2 = var_z_bzx*rc2/rzx2
#	c = normal(0, sigmac_2, n)	
	if not c:
		c = normal(0, sigmac_2)
#	ezx = normal(0, var_w_bzx*(1/rzx2 -1)-sigmac_2, n)
	ezx = normal(0, var_w_bzx*(1/rzx2 -1)-sigmac_2)
#	print(z_bzx.shape)
#	print(c.shape)
#	print(ezx.shape)
	x = z_bzx + c + ezx
	x = sum(x)
	lst = list(r1) + list(r2)
	geno_012 = list(w) + list(w2)
	print(len(geno_012))
	print(len(lst))
#	causal = dict(zip(lst, geno_012))
	causal = pd.DataFrame()
	causal['012'] = geno_012
	causal['SNP'] = lst
	ref='geno_ref.txt'
	ref = pd.read_csv(ref)
	geno_GCTA = pd.merge(causal, ref, left_on='SNP', right_on='SNP')
#	print(geno_GCTA.head())
	pool = Pool(30)
	r = []
	r2 = []
	for i in range(geno_GCTA.shape[0]):	
		r.append(pool.apply_async(gene_expr_worker, [geno_GCTA.iloc[i, :]],))
	for i in r:
		r2.append(i.get())
	geno_GCTA = pd.concat(r2, axis=1).T
	print(geno_GCTA.head())
#	geno_GCTA = geno_GCTA.apply(gene_expr_worker, axis=1)
	geno_GCTA = geno_GCTA.sort_values(by='SNP')
	print(geno_GCTA.head())
	out = ''.join(geno_GCTA['new'])
	print(x, out[:10])
	return (x, out, sigmac_2, c)

# direct causality model	
def simu_3(n, m, ind, causal=False, causal_ind=0, test=0): 
	print('into simu_3')
	rzx2 = 0.05
	rxy2 = 5e-03
	bzx = 0.6
	bxy = 0.1
	args={'rzx2':rzx2, 'rxy2':rxy2, 'bzx': bzx, 'bxy':bxy}
	# n: number of causal snps	
	# m: number of all snps
	# number of causal gene
	gene_num = np.random.randint(5, 10)
	gene_num = 2
	print(gene_num)
	x = []
	c_lst = []
	sigmac_2_lst = []
	c = 0 
	count = 0
	for i in range(gene_num):
		print('i is ', i)
		if test == 0:
			print('test 0')
			tmp_x, geno, sigmac_2, c = gene_expr(causal_lst=True, args=args)
			break
		elif test == 1:
			print('test 1')
			tmp_x, geno, sigmac_2, c = gene_expr(causal_files=['tmp/r1_0', 'tmp/r2_0', 'tmp/p_0', 'tmp/w_0', 'tmp/w2_0'], args=args)
			break
		else:
			print('normal running...')
			if i == 0:
				while 1:
					try:
						tmp_x, geno, sigmac_2, tmp_c = gene_expr(causal_lst=True, args=args)		
					#	print('c is ', c)
					#	c_lst.append(c)
						c += sum(tmp_c)
						sigmac_2_lst.append(sigmac_2)
						count += 1
						x.append(tmp_x)
					except:
						pass
					if count == 1:
						break
				
			else:
				while 1:
					try:
						tmp_x, geno, sigmac_2, c = gene_expr(causal_files=['tmp/r1_0', 'tmp/r2_0', 'tmp/p_0', 'tmp/w_0', 'tmp/w2_0'], args=args)
					#	print('c is ', c)
						c_lst.append(c)
						sigmac_2_lst.append(sigmac_2)
						count += 1
						x.append(tmp_x)
					except:
						pass
					if count == gene_num:
						break
			
	print(count)			
	x = np.array(x)	
	print('x is', x)			
	x_bxy = x*bxy
	print('x_bxy shape', x_bxy.shape)
	var_x_bxy = x_bxy.var()
	exy = []
	for sigmac_2 in sigmac_2_lst:
		exy.append(normal(0, var_x_bxy*(1/rxy2-1)-sigmac_2, n))
	c = []
	for i in c_lst:
		c.append(i)
	y = x*bxy + c + exy
	print(x)
	print(y)
	return sum(y)
if __name__ == '__main__':	
	# test
	simu_3(10000, 100000, 0, causal=False, test=2)
	raise
	## repeate 10000 times
	count = 0
	while True:
		print('into main')
		try:
			with open('phe_{}'.format(count), 'w+') as f:
				y = simu_3(10000, 100000, 0, causal=False, test=False)
				print(y)
			f.write(str(y)+'\n')
			count += 1
		except:
			pass
		if count > 1:
			break	
