import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
flowering_gene_TAIR = '/home/lhjiang/project/GWAS_eQTL/arabidopsis/data/flowering_gene'
all_gene_file = '/home/lhjiang/project/GWAS_eQTL/arabidopsis/data/all_gene'
def _stat(TP,FP,TN,FN, out_file):
        FDR = FP/(TP+FP)
        FNR = FN/(TP+FN)
        FPR = FP/(FP+TN)
        TPR = TP/(TP+FN)
        TNR = TN/(TN+FP)
        Acc = (TP+TN)/(TP+FP+TN+FN)    
        print("TP, TN, FP, FN".format(TP, TN, FP, FN))
        print('true positive rate: {}'.format(TPR))
        print('false positive rate: {}'.format(FPR))
        print('true negtive rate: {}'.format(TNR))
        print('false negtive rate: {}'.format(FNR))
        print('accuracy is {}'.format(Acc))
        with open(out_file, 'w+') as f:
                f.write('TP is {}'.format(TP)+'\n')
                f.write('FP is {}'.format(FP)+'\n')
                f.write('FN is {}'.format(FN)+'\n')
                f.write('TN is {}'.format(TN)+'\n')
                f.write('true positive rate: {}'.format(TPR)+'\n')
                f.write('false positive rate: {}'.format(FPR)+'\n')
                f.write('true negtive rate: {}'.format(TNR)+'\n')
                f.write('false negtive rate: {}'.format(FNR)+'\n')
                f.write('accuracy is {}'.format(Acc))
def confusion_mat(gene_file, out_file):
        flowering_gene_table = pd.read_csv(flowering_gene_TAIR, sep='\t')
        flowering_gene = list(flowering_gene_table['locus_name'])
        gene = open(gene_file).readlines()
        gene = [i.strip() for i in gene]
        gene = list(set(gene))
        all_gene = open(all_gene_file).readlines()
        all_gene = [i.strip() for i in all_gene]
        all_gene = list(set(all_gene))
        TP = len(set(gene) & set(flowering_gene))
        FP = len(set(gene) - set(flowering_gene))
        FN = len(set(set(all_gene) - set(gene)) & set(flowering_gene))
        TN = len(set(set(all_gene) - set(flowering_gene)) & set(set(all_gene)-set(gene)))

        TP = len(set(gene) & set(flowering_gene))
        FP = len(set(gene) - set(flowering_gene))
        FN = len(set(set(all_gene) - set(gene)) & set(flowering_gene))
        TN = len(set(set(all_gene) - set(flowering_gene)) & set(set(all_gene)-set(gene)))
        _stat(TP,FP,TN,FN, out_file)
def label_backup(g_lst):
	g_tb = pd.DataFrame(g_f, header=None)
	g_tb.columns = ['gene']
	f_tb = pd.DataFrame(f_f, header=None)	
	f_tb.columns = ['gene']
	f_tb['true'] = 1
	tb = pd.merge(g_tb, f_tb, left_on='gene', right_on='gene', how='left')
	tb = tb.fillna(0)
	tmp_tb = pd.DataFrame(g_lst)
	tmp_tb.columns = ['gene']
	tmp_tb['pred'] = 1
	all_tb = pd.merge(tb, tmp, left_on='gene', right_on='gene', how='left')
	all_tb.fillna(0)
	pred = all_tb['true']
	true = all_tb['pred']
	return true, pred
def label(f, g_f='data/arabi.gene', f_f='data/flowering_gene'): # f: two cols; gene score; 
	g_tb = pd.read_csv(g_f, header=None)
	g_tb.columns = ['gene']
	f_tb = pd.read_csv(f_f, header=None)	
	f_tb.columns = ['gene']
	f_tb['true'] = 1
	tb = pd.merge(g_tb, f_tb, left_on='gene', right_on='gene', how='left')
	tb = tb.fillna(0)
	tmp_tb = pd.read_csv(f, header=None, sep=' ')
	tmp_tb.columns = ['gene', 'score']
	tmp_tb = tmp_tb.replace(np.inf, np.nan)
	tmp_tb['score'] = tmp_tb['score'].fillna(tmp_tb['score'].max())
	score =	np.expand_dims(tmp_tb['score'].values, axis=1)
	scaler = MinMaxScaler()
#	return score, scaler
	score = scaler.fit_transform(score)
	tmp_tb['score'] = score
	all_tb = pd.merge(tb, tmp_tb, left_on='gene', right_on='gene', how='left')
	all_tb = all_tb.fillna(0)
	score = all_tb['score'].values
	true = all_tb['true'].values
	return true, score
		
