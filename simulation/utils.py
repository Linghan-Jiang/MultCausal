import os
import pandas as pd
import numpy as np
from multiprocessing import Pool
# A1: minor
# A2: major
# two minor: 2
# two major: 0
# one minor, one major: 1

f = 'test.frq'
def make_ref_worker(x):
#	print(type(x))
#	print(x)
#	print(x['A1'])
	if type(x) == list:
		x = x[0]
	x['2'] = ' ' + x['A1'] + ' ' + x['A1']
	x['1'] = ' ' + x['A1'] + ' ' + x['A2']
	x['0'] = ' ' + x['A2'] + ' ' + x['A2']
	return x
def make_ref(test=False):
	# test
	if test:
		d = pd.read_csv(f, sep='\s*', nrows=3)
#		return d
		print(d)
		d = d.iloc[:5, :].apply(make_ref_worker, axis=1)	
		print(d)
	else:
		d = pd.read_csv(f, sep='\s*')
		pool = Pool(10)
		r = []
		r2 = []
		for i in range(d.shape[0]):
		#	r.append(pool.apply_async(make_ref_worker, [d.iloc[i, :]],))
			r2.append(make_ref_worker(d.iloc[i, :]))
			
		for i in r:
			r2.append(i.get())	
		d = pd.concat(r2, axis=1)
#		d.apply(make_ref_worker, axis=1)
	d.T.to_csv('geno_ref.txt', index=False)
	return d
def num_2_ped(snp_012, snp):
	snp_012 = [str(i) for i in snp_012]
	ref = pd.read_csv('geno_ref.txt')
	tmp = ref[ref['SNP'].isin(snp)][['SNP', '0', '1', '2']]
	geno = pd.DataFrame()
	r1 = []
	r2 = []
	for i in ref.shape[0]:
		r1.append(tmp.iloc[i, 'SNP'])
		r2.append(tmp.iloc[i, snp_012[i]])
		print(r1)
		print(r2)
		break
	geno['SNP'] = r1
	geno['012'] = r2
	geno = geno.sort_by('SNP')
	result = ' '.join(list(geno['012']))
	print(result)
	return result
def bim_filter():
	bim = pd.read_csv('arabi_1135.bim', sep=' ', header=None)
	frq = pd.read_csv('arabi_1135.frq', sep='\s*', header=None)
	
	bim[frq]
def sample(n, ind=0):
	f='tmp/sample_{}_{}'.format(ind,n)
	a = pd.read_csv('arabi_1135.fam', sep=' ', header=None)
	tmp = np.random.choice(a.index, n)
	ind = a.loc[tmp, 0]
	if a.shape[1] > 1:
		fam = a.loc[tmp, 1]
	with open(f, 'w+') as f:
		for i,j in zip(fam, ind):
			f.write(str(i) + ' ' + str(j) + '\n')
	return a[[0,1]].values
def snp(n, ind=0):
	f='tmp/snp_{}_{}'.format(ind,n)
	tmp = pd.read_csv('arabi_1135_chr1_100kbwindow_15mb_0.5r2.prune.out', header=None)
	snp = np.random.choice(tmp[0], n)
	with open(f, 'w+') as f:
		for i in snp:
			f.write(str(i) + '\n')
	return snp
if __name__ == '__main__':
#	r = make_ref(False)
	sample(10)
