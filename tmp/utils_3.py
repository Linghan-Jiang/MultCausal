# utils_2.py
# This script is for the programs other than my own MatMult method.
import os
import numpy as np
import pandas as pd
import subprocess
import scipy.io
from utils import sort
from multiprocessing import Pool
def dada_input(gdic_f='data/arabi_gdic.txt', ppi='data/ppi/atpin.ppi.2evidence.no-weight.txt', assoc='data/gsmr/outcome/arabi.qassoc', sg_f='data/arabi.20kb.sg', gene_f='data/DADA.candidate.txt', ppi_f='data/DADA.ppi.txt', seed_f='data/DADA.seed.txt'): # gdic_f must derive from expression file.
	# ppi_f
	sort(ppi, gdic_f, out=ppi_f)
	tb_pp = pd.read_csv(ppi_f, header=None, sep=' ')
	tb_pp.columns = ['gene_1', 'gene_2', 'score']
	tb_ps = pd.DataFrame()
	tb_ps['gene'] = list(tb_pp['gene_1']) + list(tb_pp['gene_2'])
	# seed_f
	tb_sg = pd.read_csv(sg_f, sep=' ', header=None)
	tb_sg.columns = ['snp', 'gene']
	tb_gdic = pd.read_csv(gdic_f, sep=' ', header=None)
	tb_gdic.columns = ['ind', 'gene']
	tb_sp = pd.read_csv(assoc, sep='\s*', usecols=['SNP','P'], dtype={'P':float}, engine='python')
	tb_sp.columns = ['snp', 'p']
	tb_sgi = pd.merge(tb_sg, tb_gdic, left_on='gene', right_on='gene')
	tb_sgip = pd.merge(tb_sgi, tb_sp, left_on='snp', right_on='snp')
	tmp = tb_sgip.loc[tb_sgip.groupby('gene')['p'].idxmin()]
	tmp['p'] = np.int16(-np.log(tmp['p']))
	
	n = int(subprocess.check_output('wc -l {}'.format(gdic_f), shell=True).decode('ascii').split()[0])
	p_array = np.zeros(n)
	for i in tmp['ind']:
		ind = int(i)
		val = tmp[tmp['ind']==i]['p']
		p_array[ind] = val
	np.save(seed_f, p_array)
	
	tmp = tmp[tmp['ind'].isin(tb_ps['gene'])]
	tmp[['ind', 'p']].to_csv(seed_f, header=None, index=False, sep=' ')
	# gene_f
	tmp['ind'].to_csv(gene_f, header=None, index=False, sep=' ')
dada_input()
def gsmr_cut(f1, f2, out_f):
	if os.path.exists(out_f):
		os.remove(out_f)
	order = "echo 'SNP A1 A2 freq b se p N' >> {}".format(out_f); os.system(order)
	order = "awk '(NR==FNR)&&(NR>1){{a[$2]=$2 FS $3 FS $4 FS $5;next}}($2 in a){{print a[$2], $5, $6, $9, $4}}' {} {} >> {}".format(f1, f2, out_f); os.system(order)
def gsmr_baby(g_f, h_f, out_d):
        with open(g_f, 'w+') as f:
                pool = Pool(30)
                for i in os.listdir(out_d):
                        if i.endswith('qassoc'):
                                tmp_f = os.path.abspath(os.path.join(out_d, i)) 
                                tmp_out = tmp_f.replace('qassoc', 'gsmr')
                                pool.apply_async(gsmr_cut, (h_f, tmp_f, tmp_out,))    
                                f.write('{} {}\n'.format(i.split('.')[0], tmp_out))
                pool.close()
                pool.join()
def gsmr_file(b_pre='test/test', exp_f='test/test.plink.exp', out_d='./test_gsmr', n=10, test=None):
	exposure_d = '{}/exposure'.format(out_d)
	outcome_d = '{}/outcome'.format(out_d)
	if not os.path.isdir(exposure_d):
		os.makedirs(exposure_d)
	if not os.path.isdir(outcome_d):
		os.makedirs(outcome_d)
	exposure_f = '{}/gsmr_exposure.txt'.format(out_d)
	outcome_f = '{}/gsmr_outcome.txt'.format(out_d)
	g_lst = open(exp_f).readline().strip().split()[2:]
	if test:
		g_lst = g_lst[:1]
		n = 1 
	# exposure
	pool = Pool(n)
	for g in g_lst:
		order = 'plink --bfile {} --pheno {} --geno 0.05 --pheno-name {} --assoc --allow-no-sex --out {}/{}'.format(b_pre, exp_f, g, exposure_d, g); print(order)
		pool.apply_async(os.system, (order,))
	pool.close()
	pool.join()
	# outcome
	order = 'plink --bfile {} --geno 0.05 --assoc --allow-no-sex --out {}/arabi'.format(b_pre, outcome_d); os.system(order)
	# h_f
	order = 'plink --bfile {} --freq --out {}/arabi'.format(b_pre, out_d); os.system(order)
	h_f = '{}/arabi.frq'.format(out_d)
	# gsmr cut
	gsmr_baby(exposure_f, h_f, exposure_d)
	gsmr_baby(outcome_f, h_f, outcome_d)
def chr():
        tb['chr'][8428147:10709949] = 5
def dada_out(f='/home/lhjiang/software/network_tools/dada/source/matlab.mat', out_f=None):
	import scipy.io
	if not out_f:
		out_f = f[:-4]+'.csv'
	tmp = scipy.io.loadmat(f)
	tmp = pd.DataFrame(tmp['proteinList'].transpose())
	tmp.apply(lambda x: x[0][0], axis=1).to_csv(out_f, header=None, index=False)
def dada_input(gdic_f='data/gdic.txt', ppi='data/atpin.ppi.2evidence.no-weight.txt', assoc='results/outcome/test.qassoc', sg_f='data/arabi.20kb.sg', gene_f='data/DADA.candidate.txt', ppi_f='data/DADA.ppi.txt', seed_f='data/DADA.seed.txt'): # gdic_f must derive from expression file.
	# ppi_f
	sort(ppi, gdic_f, out=ppi_f)
	tb_pp = pd.read_csv(ppi_f, header=None, sep=' ')
	tb_pp.columns = ['gene_1', 'gene_2', 'score']
	tb_ps = pd.DataFrame()
	tb_ps['gene'] = list(tb_pp['gene_1']) + list(tb_pp['gene_2'])
	# seed_f
	tb_sg = pd.read_csv(sg_f, sep=' ', header=None)
	tb_sg.columns = ['snp', 'gene']
	tb_gdic = pd.read_csv(gdic_f, sep=' ', header=None)
	tb_gdic.columns = ['ind', 'gene']
	tb_sp = pd.read_csv(assoc, sep='\s*', usecols=['SNP','P'], dtype={'P':float}, engine='python')
	tb_sp.columns = ['snp', 'p']
	tb_sgi = pd.merge(tb_sg, tb_gdic, left_on='gene', right_on='gene')
	tb_sgip = pd.merge(tb_sgi, tb_sp, left_on='snp', right_on='snp')
	tmp = tb_sgip.loc[tb_sgip.groupby('gene')['p'].idxmin()]
	tmp['p'] = np.int16(-np.log(tmp['p']))
	
	n = int(subprocess.check_output('wc -l {}'.format(gdic_f), shell=True).decode('ascii').split()[0])
	p_array = np.zeros(n)
	for i in tmp['ind']:
		ind = int(i)
		val = tmp[tmp['ind']==i]['p']
		p_array[ind] = val
	np.save(seed_f, p_array)
	
	tmp = tmp[tmp['ind'].isin(tb_ps['gene'])]
	tmp[['ind', 'p']].to_csv(seed_f, header=None, index=False, sep=' ')
	# gene_f
	tmp['ind'].to_csv(gene_f, header=None, index=False, sep=' ')
def simu_phe(fam='./test/gsmr_example.fam', out_pre='./test/gsmr_example'):
	tmp = pd.read_csv(fam, sep=' ', header=None)
	tmp.columns = ['FID', 'IID', 'FATHER', 'MATHER', 'SEX', 'PHE']
	tmp['PHE'] = np.random.normal(0, 0.1, tmp.shape[0])
	tmp['C1'] = np.random.normal(0, 0.1, tmp.shape[0])
	tmp['C2'] = np.random.normal(0, 0.1, tmp.shape[0])
	
	out_f1 = out_pre+'.phe'
	out_f2 = out_pre+'.exp'
	tmp[['FID', 'IID', 'PHE']].to_csv(out_f1, index=False, sep=' ')
	tmp[['FID', 'IID', 'C1', 'C2']].to_csv(out_f2, index=False, sep=' ')
def corr(pl_f='./test_2/plink.qassoc', ge_f='./test_2/gemma.assoc.txt'):
	pl = pd.read_csv(pl_f, sep='\s*')
	ge = pd.read_csv(ge_f, sep='\s*')
	tmp = pd.merge(pl, ge, left_on='SNP', right_on='rs')
	tmp = tmp[['P', 'p_wald']].dropna().values
	print(tmp.shape)
	print(np.corrcoef(tmp[:,0], tmp[:,1]))
def plink_exp(fam='./test/test.fam', exp_f='./test/test_2.exp', out_d=None, g_f=None): #
	if not out_d:
		out_d, pre = os.path.split(exp_f)
	else:
		tmp_d, pre = os.path.split(exp_f)
	fam = pd.read_csv(fam, sep=' ', header=None, dtype={1:str})
	fam.columns = ['FID', 'IID', 'FATHER', 'MATHER', 'SEX', 'PHE']
	fam.index = fam['IID']
	exp = pd.read_csv(exp_f, sep=' ', index_col=0).T
	if g_f:	
		cols = ['FID', 'IID'] + [i.strip() for i in open(g_f).readlines()]
	else:
		cols = ['FID', 'IID'] + list(exp.columns)
	print(fam.iloc[:5, :5])
	print(exp.iloc[:5, :5])
	print(fam.index.dtype)
	print(exp.index.dtype)
	tmp = pd.merge(exp, fam, left_index=True, right_index=True, how='left')
	print(tmp.iloc[:5, :5])
	pre = pre.split('.')[0]
	out_f = '{}/{}.plink.exp'.format(out_d, pre)
	tmp['FID'] = tmp.index
	tmp['IID'] = tmp.index
	tmp[cols].to_csv(out_f, na_rep='NA', sep=' ', index=False)
	return exp
#plink_exp(fam='data/arabi.fam', exp_f='data/GSE54680_10C_quant.csv', out_d='data', g_f='data/arabi.gene')
def gemma_input(b_pre, phe=None, out_pre=None):
	if not out_pre:
		out_pre = b_pre
		
	fam_f = b_pre + '.fam'
	pheno_f = out_pre + '.gemma.pheno'
	geno_f = out_pre + '.gemma.geno'
	annot_f = out_pre + '.gemma.annot'
	tmp_0 = out_pre + '.gemma.traw'
	# pheno
	if not phe:
		order = "awk '{{print $6}}' {} > {}".format(fam_f, pheno_f); os.system(order)
	else:
		order = """awk '{{out=""; for(i=3;i<=NF;i++){{out=out" "$i}}; print out}}' {} > {}""".format(phe, pheno_f); os.system(order)
	# snp
	order = "plink --bfile {} --recode A-transpose --out {}.gemma".format(b_pre, out_pre); os.system(order)
	order = """awk '(NR>1){{OFS=","; out=$2 OFS $5 OFS $6; for(i=7;i<=NF;i++){{out=out OFS $i}}; print out}}' {} > {}""".format(tmp_0, geno_f); os.system(order)
	# annot 
	order = "awk '(NR>1){{print $2, $4, $1}}' {} > {}".format(tmp_0, annot_f); os.system(order)
	# relationship
	if '/' in out_pre:
		out_dir, out_pre = os.path.split(out_pre)
	else:
		out_dir = './'
	order = 'gemma -g {} -p {} -o {} -outdir {} -gk 1'.format(geno_f, pheno_f, out_pre, out_dir); print(order);os.system(order)
def matrixeqtl_input(b_pre='test/test', exp='data/GSE54680_10C_raw_quantile.txt', out_pre='test/test_2'):
	if not out_pre:
		out_pre = b_pre
	fam_f = b_pre + '.fam'
	exp_f = out_pre + '.exp'
	snp_f = out_pre + '.snp'
	snploc_f = out_pre + '.snploc'
	geneloc_f = out_pre + '.geneloc'
	keep = out_pre + 'keep'
	snp_tmp = out_pre + '.traw'
	# common
	col_1 = open(exp).readline().strip().split()
	col_2 = list(pd.read_csv(fam_f, sep=' ', usecols=[0], header=None, dtype={0:str}).values.flatten())
	cols = sorted([int(i) for i in (set(col_1) & set(col_2))])
	print('length of indivisual ', len(cols))
	with open(keep, 'w+') as f:
		for c in cols:
			f.write('{} {}\n'.format(c, c))
	# exp
	exp = pd.read_csv(exp, sep=' ', index_col=0)
	exp.columns = [int(i) for i in exp.columns]
	exp = exp.dropna(how='all')
	if os.path.exists(exp_f):
		print('WARNING: {} ALLREADY EXISTED.'.format(exp_f)) 
	exp[cols].to_csv(exp_f, sep=' ', na_rep='NA')
	# snp
	print('MAKING SNP FILE')
	order = "plink --bfile {} --recode A-transpose --keep {} --out {}".format(b_pre, keep, out_pre); os.system(order)
	order = """awk '{{out=$2; for(i=7;i<=NF;i++){{out=out" "$i}}; print out}}' {} > {}""".format(snp_tmp, snp_f); os.system(order)
	header = [int(i.split('_')[0]) for i in open(snp_tmp).readline().strip().split()[6:]]
	if header != cols:
		raise('INDIVISUAL ORDER ERROR!')
	header = 'snpid ' + ' '.join([str(i) for i in header])
	order = "sed -i '1s/^.*/{}/g' {}".format(header, snp_f); os.system(order)
	print('SNP FILE DONE')
	# snploc
	print('MAKING SNPLOC FILE')
	order = "awk '(NR>1){{print $2, $1, $4}}' {} > {}".format(snp_tmp, snploc_f); os.system(order)
	print('SNPLOC FILE DONE')
	order = 'rm {}.traw {}.log {}.nosex'.format(out_pre, out_pre, out_pre); os.system(order)
if __name__ == '__main__':
	#gsmr_cut('test/gsmr_2/test.frq', 'test/gsmr_2/test.qassoc', 'test/gsmr_2/test.phe.gsmr')
	#gemma_input('data/arabi', out_pre='data/gemma/arabi')	
	pass	
