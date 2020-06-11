import os
import sys
import scipy.io
import scipy.sparse
from utils_2 import *
from multiprocessing import Pool
def gsmr(bfile='test_2/test', out='results/gsmr_5e-08/test', exposure='./results/exposure', outcome='./results/outcome', gwas_thresh=5e-08):
	order = 'gcta64 --bfile {} --gsmr-file {} {} --gsmr-direction 0 --out {}'.format(b_pre, f1, f2, out_pre)
	os.system(order)
def gene_rank(ppi_f='data/DADA.ppi.txt.npz', p_f='data/DADA.seed.txt.npy', d=0.3,test=None, out='results/gene_rank_0.3.txt', gdic_f='data/arabi_gdic.txt'):
	if test:
		tmp = scipy.io.loadmat(test)
		w = tmp['w_All']
		p = tmp['expr_data']
	else:
		w = scipy.sparse.load_npz(ppi_f)
		p = np.load(p_f)
	p = np.abs(p)
	norm_p = p/np.amax(p, axis=0)
	degrees = np.sum(w, axis=0)
	degrees[degrees == 0] = 1
	degrees = np.array(degrees)
	d1 = scipy.sparse.csr_matrix(np.diag((1/degrees).transpose().flatten()))
	print('d1.shape: ', d1.shape)
	print('w.shape: ', w.shape)
	print('norm_p', norm_p.shape)
	print('type w', type(w))
	print('type d1', type(d1))
	print('type d', type(d))
	a = scipy.sparse.csc_matrix(np.eye(w.shape[0]) - d * (w.transpose().dot(d1)))
	print('type a', type(a))

	b = (1-d) * norm_p
	r = scipy.sparse.linalg.inv(a).dot(b)
	r_tb = pd.DataFrame(r)
	r_tb.columns = ['p']
	gdic_tb = pd.read_csv(gdic_f, sep=' ', header=None, index_col=0)
	gdic_tb.columns = ['gene']
	tmp = pd.merge(r_tb, gdic_tb, left_index=True, right_index=True)
	tmp[['gene', 'p']].to_csv(out, header=None, sep=' ', index=False)
	return r
if __name__ == '__main__':
	for d in range(1, 10):
		d = 0.1 * d
		gene_rank(out='results/gene_rank/gene_rank.{}.txt'.format(d), d=d)
