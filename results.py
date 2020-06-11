# Get final results
import os
import random
from stat_jiang import label
from sklearn.metrics import roc_auc_score, roc_curve
from matplotlib import pyplot as plt
def gsmr_worker(d, out):
	tmp_f = '{}/tmp_gsmr{}'.format(d, str(random.random())[-10:])
	order = 'cat {}/*.gsmr > {}'.format(d, tmp_f)
	os.system(order)
#	print(order)
	order = """awk '(NR>1){{if ($4!="nan" && $4!="-nan")print $1, -log($5)}}' {} > {}""".format(tmp_f, out)
#	print(order)
	os.system(order)
	os.remove(tmp_f)
def gemma(f, out):
	order = "awk '(NR>1){{print $2, -log($11)}}' {} > {}".format(f, out)
#	print(order)
	os.system(order)
def gsmr():
 	## gsmr
	for i in range(2, 9):
		gsmr_worker('data/gsmr/gsmr_5e-0{}'.format(i), 'results/gsmr/gsmr_5e-0{}.txt'.format(i))
	raise
def mult():
	## test
#	gemma('data/mult/arabi.0.01.1e+06Dist.p.atpin.2.weight.ppi.isg.assoc.txt', 'results/mult/arabi.0.01.1e+06Dist.p.atpin.2.weight.ppi.isg.txt')
	## mult
	mult_d = 'data/mult'
	out_d = 'results/mult'
	for i in os.listdir(mult_d):
		if i.endswith('assoc.txt') and i.startswith('arabi.20kb.test'):
			f = '{}/{}'.format(mult_d, i)
			out = '{}/{}.txt'.format(out_d, i[:-10])
#			print(f, out)
			gemma(f, out)
def main(d, fig=None, method=None):	
	dic = dict()
	for i in os.listdir(d):
		if i.endswith('txt'):#and ('20kb' in i):
			f = '{}/{}'.format(d, i)
			true, score = label(f)
			auc = roc_auc_score(true, score)
			if fig:
				fpr, tpr, thresh = roc_curve(true, score)
				dic[auc] = [fpr, tpr]
			print(i.rjust(60), auc)
	tmp = sorted(dic)[-1]
	fpr = dic[tmp][0]
	tpr = dic[tmp][1]
	plt.plot(fpr, fpr, color='grey')
	plt.plot(fpr, tpr, label='{}, auc={:.2f}'.format(method, tmp))
	return fig

if __name__ == '__main__':
	d_1 = 'results/gsmr'
	d_2 = 'results/mult'
	d_3 = 'results/gene_rank'
	d_4 = 'results/gwas'
#	mult()
#	main(d_1)
#	main(d_4)
#	main(d_3)
	fig = plt.figure()
	fig = main(d_2, fig, 'MultCausal')
	fig = main(d_3, fig, 'GeneRank')
	fig = main(d_1, fig, 'GSMR')
	fig = main(d_4, fig, 'LM')
	plt.legend()
	plt.xlabel('False Positive Rate')	
	plt.ylabel('True Positive Rate')	
#	plt.show()
	plt.savefig('test.png')
	
