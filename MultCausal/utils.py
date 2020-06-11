import os
import scipy
import random
import subprocess
import numpy as np
import pandas as pd
from multiprocessing import Pool
from scipy.sparse import csr_matrix, lil_matrix
from sklearn.preprocessing import MinMaxScaler
def bimbam_2_is_baby(f, out_f):
	simat = pd.read_csv(f, header=None)
	simat = simat.iloc[:, 3:].values
	ismat = csr_matrix(simat).transpose()	
	scipy.sparse.save_npz(out_f, ismat)
def bimbam_2_is(f, si_d, is_d, step=1000):
	print('Into bimbam_2_is')
	order = 'split -d -a 5 -l {} {} {}/'.format(step, f, si_d); os.system(order)
	pool = Pool(10)
	for i in os.listdir(si_d):
		f = os.path.join(si_d, i)
		out_f = os.path.join(is_d, i)
		pool.apply_async(bimbam_2_is_baby, (f, out_f,))
	pool.close()
	pool.join()	
def mateqtl_2_sg_baby(f, method):
	d, f_tmp = os.path.split(f)
	d = os.path.join(d, method)
	out_f = os.path.join(d, f_tmp)
	if method == 'p':
		order = "awk '{{print $1, $2, -log($5)}}' {} > {}".format(f, out_f); os.system(order)
	elif method == 'beta':
		order = "awk '{{print $1, $2, $3}}' {} > {}".format(f, out_f); os.system(order)
	elif method == 'abs-beta':
		order = "awk '{{if ($3>0) a=$3; else a=-$3;print $1, $2, a}}' {} > {}".format(f, out_f); os.system(order)
	elif method == 'both':
		order = "awk '{{print $1, $2, -log($5)*$3}}' {} > {}".format(f, out_f); os.system(order)
	elif method =='abs-both':
		order = "awk '{{if ($3>0) a=$3; else a=-$3;print $1, $2, -log($5)*a}}' {} > {}".format(f, out_f); os.system(order)		
def mateqtl_prepare(pre):
	cis = pre + '.cis'
	tra = pre + '.tra'
	tmp_d = 'tmp_' + str(random.random())[-10:]
	if not os.path.isdir(tmp_d):
		os.makedirs(tmp_d)
	if ('SNP' in open(cis).readline()) and ('gene' in open(cis).readline()):		
		order = "sed -i '1d' {}".format(cis); os.system(order)
	if ('SNP' in open(tra).readline()) and ('gene' in open(tra).readline()):		
		order = "sed -i '1d' {}".format(tra); os.system(order)
	order = 'split -l 10000 -a 5 {} {}/'.format(tra, tmp_d); os.system(order); print(order)
	order = 'split -l 10000 -a 4 {} {}/'.format(cis, tmp_d); os.system(order)
	print('Prepare done'.ljust(60, '-'))
	return tmp_d
def mateqtl_2_sg(pre='./test/test', method='p', tmp_d=None): # Alternertives: 'both', 'beta', 'abs-beta', 'abs-both'
	print('Perfoming mateqtl_2_sg.')
	if not tmp_d:
		tmp_d = mateqtl_prepare(pre)
	print('tmp_d is {}'.format(tmp_d).ljust(60, '-'))
	out_d = os.path.join(tmp_d, method)
	if not os.path.isdir(out_d):
		os.makedirs(out_d)
	pool = Pool(20)
	for i in os.listdir(tmp_d):
		i = os.path.join(tmp_d, i)
		if os.path.isfile(i):
			pool.apply_async(mateqtl_2_sg_baby, (i, method,))
	pool.close()
	pool.join()
	d, tmp_pre = os.path.split(pre)
	out_f = '{}/{}.{}.sg'.format(d, tmp_pre, method)
	order = 'cat {}/{}/* > {}'.format(tmp_d, method, out_f); os.system(order);# print(order)
	print('mateqtl_2_sg for {} done !'.format(pre).ljust(60, '-'))
	print('result saved to {}'.format(out_f).ljust(60, '-'))
	return tmp_d
def gsmr_exp_sort(sg_f, out_d, ind):
	tmp = pd.read_csv(sg_f, sep='\t')
	tmp.columns = ['SNP', 'gene', 'beta', 't-stat', 'p-value', 'FDR']
	tmp = tmp.replace(np.inf, np.nan)
	tmp['t-stat'] = tmp['t-stat'].fillna(tmp['t-stat'].max())
	tmp['se'] = np.sqrt(tmp['beta']/tmp['t-stat'])	
	grps = tmp.groupby('gene')
	for g in grps:
		out_f = '{}_{}'.format(g[0], ind)
		out_f =	os.path.join(out_d, out_f)
		out(g[1], out_f)
def gsmr_exp_merge(g, d, tmp_d, o_f): # d is gsmr_d
	f = os.path.join(d, g)
	for i in os.listdir(tmp_d):
		if i.startswith(g):
			order = 'awk (NR==FNR){{a[$1]}}($2 in a){{print }}'
			order = "cat {} >> {} ".format(os.path.join(tmp_d, i), f); os.system(order) 
def gsmr_exp(pre='./test/test', gsmr_d='./test/gsmr', tmp_d=None, n=20, g_f='./test/test_2.gene'):
	gsmr_tmp = os.path.join(gsmr_d, 'tmp')
	print('tmp dir is {}'.format(gsmr_tmp).ljust(60, '-'))
	# sort_by_gene
	if not tmp_d:
		tmp_d = mateqtl_prepare(pre)
	if not os.path.isdir(gsmr_d):
		os.makedirs(gsmr_d)
	if not os.path.isdir(gsmr_tmp):
		os.makedirs(gsmr_tmp)
	print('Working using {} cores'.format(n).ljust(60, '-'))
	print('Sorting begins'.ljust(60, '-'))
	pool = Pool(n)
	for ind, f in enumerate(os.listdir(tmp_d)):
		f = os.path.join(tmp_d, f)
		gsmr_exp_sort(f, gsmr_tmp, ind)
#		raise
#		break
#		pool.apply_async(gsmr_exp_sort, (f, tmp_d, ind,))
	pool.close()
	pool.join()
	print('Sorting done'.ljust(60, '-'))
	print('Merging begins'.ljust(60, '-'))
	g_lst = [i.strip() for i in open(g_f).readlines()]
	pool = Pool(n)
	for g in g_lst:
		pool.apply_async(gsmr_exp_merge, (g, gsmr_d, tmp_d,))
	pool.close()
	pool.join()
	print('Merging done'.ljust(60, '-'))
	print('result saved to {}'.format(gsmr_d).ljust(60, '-'))
def left(gdic_f):
	tmp = pd.read_csv(gdic_f, sep=' ', header=None)
	frame = pd.DataFrame()
	frame['0'] = tmp.iloc[:,1]
	frame['1'] = [random.choice(['A','T','C','G']) for i in range(tmp.shape[0])]
	frame['2'] = [random.choice(['A','T','C','G']) for i in range(tmp.shape[0])]
	return frame
def rescale(mat):
        scaler = MinMaxScaler()
        scaler.set_params(**{'feature_range':(0, 2)})
        mat = scaler.fit_transform(mat)
        return mat 
def sort(mat_f, dic_f, dic_f2=None, out=None, dtype='int', rm=None):
	# mat into number format
	if out:
		out_tmp = out
		out_npz = '{}.npz'.format(out)
		if os.path.exists(out):
			os.remove(out)
			print('The out file exists, it will be overlaped.')
	else:
		out_tmp = str(random.random())[-10:]
	if dic_f2:
		order = """awk '(FILENAME==ARGV[1]){{a[$2]=$1}}(FILENAME==ARGV[2]){{b[$2]=$1}}(FILENAME==ARGV[3])&&($1 in a)&&($2 in b){{print a[$1], b[$2], $3}}' {} {} {} > {}""".format(dic_f, dic_f2, mat_f, out_tmp)
		n = int(subprocess.check_output('wc -l {}'.format(dic_f2), shell=True).split()[0])
	else:
		order = """awk '(NR==FNR){{a[$2]=$1}}($1 in a)&&($2 in a){{print a[$1], a[$2], $3}}' {} {} > {} """.format(dic_f, mat_f, out_tmp)
		order = """awk '(NR==FNR){{a[$2]=$1}}($1 in a)&&($2 in a){{print a[$2], a[$1], $3}}' {} {} >> {} """.format(dic_f, mat_f, out_tmp)
		n = None
	os.system(order)
	m = int(subprocess.check_output('wc -l {}'.format(dic_f), shell=True).split()[0])
	data = pd.read_csv(out_tmp, sep=' ', header=None)
	if rm:
		os.system('rm {}'.format(out_tmp))
	if m and n:
		mat = lil_matrix((m, n), dtype=dtype)
	else:
		mat = lil_matrix((m, m), dtype=dtype)
	for i, j, k in data.values:
		mat[i, j] = k
	if out:
		mat = csr_matrix(mat)
		scipy.sparse.save_npz(out_npz, mat)
	return mat
def snp_range(bim, out_dir):
	out = str(random.random())[-10:]
	order = """awk '{{a = $4-20000; if (a > 0) print $2, $1, a, $4+20000;\
		else print $2, $1, 0, $4+20000}}' {} > {} """.format(bim, out)
	os.system(order)
	order = 'split -d -a 3 -l 1000 {} {}'.format(out, out_dir)
	os.system(order)
def gene_range(g_f, out_f=None, kb=20):
	if not out_f:
		out_f = '{}_{}kb.txt'.format(g_f.split('.')[0], kb)
	order = """awk '($3=="gene"){{b=$5+{}; c=$4-{}; match($1, "[0-9]{{1}}", d); match($9, "AT.{{7}}", e); if (d[0] != "") print e[0], d[0], c, b}}' {} > {}""".format(kb * 1000, kb *1000, g_f, out_f)
	os.system(order)
	return out_f
def sg_1_baby(bim_f, g_rng_f, out_f):
	snp = pd.read_csv(bim_f, usecols=[0, 1, 3], header=None, sep='\t', dtype={0:str, 1:str, 3:int})
	gene = pd.read_csv(g_rng_f, header=None, sep=' ', dtype={0:str, 1:str, 2:int, 3:int})
	tmp = pd.merge(snp, gene, how='left', left_on=0, right_on=1)
	tmp = tmp[tmp['3_x'] > tmp[2]]
	tmp = tmp[tmp['3_x'] < tmp['3_y']]
	tmp = tmp.dropna()
	tmp[['1_x', '0_y']].to_csv(out_f, header=None, index=False, sep=' ')
	print(tmp.iloc[:5, :])
def sg_1(bim_f, gff_f, sg_f='test/test.sg', sg_d='test/test_sg_d1'):
	if os.path.isdir(sg_d):
		raise('sg directory already existed, please choose another to place the result.')
	os.makedirs(sg_d)
	order = "split -l 10000 -a 6 {} {}/".format(bim_f, sg_d); os.system(order)
	tmp_d, _ = os.path.split(sg_d)
	g_rng_f = gene_range(gff_f)
	pool = Pool(10)
	for i in os.listdir(sg_d):
		tmp_f = '{}/{}.sg'.format(sg_d, i)
		pool.apply_async(sg_1_baby, (os.path.join(sg_d, i), g_rng_f, tmp_f, ))
	pool.close()
	pool.join()
	if os.path.exists(sg_f):
		os.system('rm {}'.format(sg_f))
	for i in os.listdir(sg_d):
		if i.endswith('sg'):
			os.system('cat {}/{} >> {}'.format(sg_d, i, sg_f))
def sg_baby(order, sg_tmp, sdic_f, gdic_f, out_f):
	os.system(order)	
	sort(sg_tmp, sdic_f, gdic_f, out=out_f, rm=True)
def sg(sg_f, sdic_f, gdic_f, si_d, sg_d):
	tmp_d = '{}_{}'.format(sg_d, str(random.random())[-10:])
	if not os.path.isdir(tmp_d):
		os.makedirs(tmp_d)
	pool = Pool(15)
	for i in os.listdir(si_d):
		si_f = '{}/{}'.format(si_d, i)
		sg_tmp = '{}/{}'.format(tmp_d, i)
		out_f = '{}/{}'.format(sg_d, i)
		order = "awk -F '[, ]' '(NR==FNR){{a[$1];next}}($1 in a){{print $0}}' {} {} > {}".format(si_f, sg_f, sg_tmp)
		pool.apply_async(sg_baby, (order, sg_tmp, sdic_f, gdic_f, out_f,))
	pool.close()
	pool.join()
def sg_2(p, d):
        tmp_d = mateqtl_2_sg('data/mateqtl/arabi.5e-0{}.1e+0{}Dist'.format(p, d), 'p')
        for m in ['both', 'beta', 'abs-beta', 'abs-both']:
                mateqtl_2_sg('data/mateqtl/arabi.5e-0{}.1e+0{}Dist'.format(p, d), m, tmp_d)
        os.system('rm -rf {}'.format(tmp_d))

if __name__ == '__main__':
#	sg_1(bim_f='test_2/test.bim', gff_f='data/TAIR10_GFF3_genes_transposons.gff', sg_f='data/arabi.20kb.sg', sg_d='data/sg_20kb')
#	for p in range(5, 9):
		p = 9
		for d in [5,6,7]:
			sg_2(p, d)
