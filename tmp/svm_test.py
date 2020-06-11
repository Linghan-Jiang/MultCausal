import os
import numpy as np
import pandas as pd
from scipy.sparse import *
from scipy.stats.stats import pearsonr
from sklearn.svm import SVR
from sklearn.metrics.pairwise import *
# sample snps
def sample_geno(f, lst, pre='data/svm/test'):
	# save snps
	with open(pre+'.snps', 'w+') as tmp:
		for i in lst:
			tmp.write(str(i)+'\n')
	tb = pd.read_csv(f)
	tb.iloc[lst, :].to_csv(pre+'.sample{}.geno'.format(len(lst)), index=False)
	x = tb.iloc[lst, 3:].T.values
	np.save(pre+'.sample{}.geno'.format(len(lst)), x)
	print(x[:5, :5])
	return x
def test_geno(size):
	lst = np.random.choice(194185, size=size)
	rows = [int(i.strip()) for i in open('data/svm/test.ind').readlines()]
	x = sample_geno('data/gemma/arabi.gemma.geno', lst)
	x = x[rows, :] 
	print(x.shape)
	return x
def exp(f='data/GSE54680_10C_quant.csv', gdic_f='data/arabi_gdic.txt', out='data/svm/test.exp', fam='data/arabi.fam'):
	iid = pd.read_csv(fam, header=None, sep=' ', usecols=[0]).values.flatten()
	iid = [str(i) for i in iid]
	lst = pd.read_csv(gdic_f, header=None, index_col=0, sep=' ').values.flatten()
	lst = [i.strip() for i in lst]
	x = pd.read_csv(f, sep=' ', index_col=0, dtype={0:str}).T	
	iid = [str(i) for i in sorted([int(i) for i in list(set(iid)&set(x.index))])]
	print(iid)
	x = x.loc[iid, lst]
	print(x.shape)
	x.to_csv(out, sep=' ')
# linear kernel
def test_1(x, pre='test', d='data/svm'):
	r = linear_kernel(x) 
	f = '{}/{}.linear.kernel'.format(d, pre)
	np.save(f, r)	
	#np.save(f, csr_matrix(f))
# polynomial kernel
def test_2(x, d='data/svm', pre='test'):
	for degree in range(2, 6):
		r = polynomial_kernel(x, degree=degree)
		f = '{}/{}.polynomial-kernel.degree{}'.format(d, pre, degree)
		np.save(f, r)
# rbf kernel
def test_3(x, d='data/svm', pre='test'):
	n_features = x.shape[1]
	gamma_lst = []
	for i in [0.1, 0.5, 1.5, 2.0, 2.5]:
		gamma_lst.append(1.0 / n_features * i)
	print('gamma_lst is : ', gamma_lst)
	for gamma in gamma_lst:
		r = rbf_kernel(x, gamma=gamma)
		f = '{}/{}.rbf-kernel.gamma{:.6f}'.format(d, pre, gamma)
		np.save(f, r)
# test above
def test():
	x = np.reshape(np.arange(6), (2,3))
	test_2(x)
	test_3(x)
# test SVR
def basic_test():
	n_samples, n_features = 10, 5
	np.random.seed(0)
	y = np.random.randn(n_samples)
	X = np.random.randn(n_samples, n_features)
	clf = SVR(gamma=1 / (n_features * X.std()), C=1.0, epsilon=0.2)
	clf.fit(X, y) 
	return r
def main(f='data/svm/test.exp', pre='test', x=None):
	if not x.any():
		print('Reading x from {}'.format(f))
		x = pd.read_csv(f, index_col=0, sep=' ').values
	print(x.shape)
	test_1(x, pre='{}.{}.{}'.format(pre, x.shape[0], x.shape[1]))
	test_2(x, pre='{}.{}.{}'.format(pre, x.shape[0], x.shape[1]))
	test_3(x, pre='{}.{}.{}'.format(pre, x.shape[0], x.shape[1]))
def ig(ig_d, idic_f='data/arab_idic.txt'):
	ind = [int(i.strip()) for i in open(idic_f).readlines()]
	ig = os.listdir(ig_d)
	igmat = scipy.sparse.load_npz(os.path.join(ig_d, ig[0]))
	for i in ig[1:]:
		igmat += scipy.sparse.load_npz(os.path.join(ig_d, i))
	igmat[ind, :].to_csv()
	return	
def compare_1(d='data/svm', ppi_f='data/ppi_mult/atpin.2.no-weight.ppi.txt.npz'):
	mat_1 = load_npz(ppi_f).toarray().flatten()
	for i in os.listdir(d):
		if i.endswith('npy'):
			mat_tmp = np.load('{}/{}'.format(d, i)).flatten()
			print(i)
			print(np.corrcoef(mat_1, mat_tmp))
	return	
def compare_2(lst_1, lst_2):
	for f1 in lst_1:
		i = np.load(f1).flatten()
		for f2 in lst_2:
			j = np.load(f2).flatten()
			print(f1, f2)
			print(pearsonr(i, j))
			
def npy_2_npz(npy_f):
	npz_f = npy_f.replace('npy', 'npz')
	tb = np.load(npy_f)	
	tb = csr_matrix(tb)
	save_npz(npz_f, tb)
def worker_1():
	for i in os.listdir('data/svm'):
		if i.endswith('npy'):
			npy_2_npz('data/svm/{}'.format(i))
	
if __name__ == '__main__':
	worker_1()
	raise
	d = 'data/svm'
	lst_0 = []
	lst_1 = []
	lst_2 = []
	for i in os.listdir(d):
		if i.startswith('test.153.1000'):
			lst_1.append('{}/{}'.format(d, i))
		elif i.startswith('test.153.2000'):
			lst_2.append('{}/{}'.format(d, i))
		elif i.startswith('test.153.23445'):
			lst_0.append('{}/{}'.format(d, i))
	compare_2(lst_1, lst_2)
	compare_2(lst_0, lst_1)
	compare_2(lst_0, lst_2)
	raise
	x = test_geno(2000)
	main(x=x)

