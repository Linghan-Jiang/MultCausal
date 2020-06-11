# Clean ppi data
import pandas as pd
import os
import networkx as nx
from utils import sort
ppi_f1 = './data/atpin_ppi.txt'
ppi_f2 = './data/string_ppi.txt'
def filt_atpin(f):
	tb = pd.read_csv(f, sep='\t')
	tb.columns = ['n1', 'n2', 'm']
	tb = tb.dropna()
	print(tb.iloc[:5, :5])
	tb['nu'] = tb['m'].apply(lambda x: x.count('|'))
	tb['w'] = 1
	tb = tb.drop_duplicates(subset=['n1','n2'])
	print(tb.shape)
	tb_2 = tb.copy(deep=True)
	n1 = list(tb['n1'])
	n2 = list(tb['n2'])
	tb_2['n1'] = n2
	tb_2['n2'] = n1
	tb = pd.concat([tb, tb_2])
	tb = tb.drop_duplicates(subset=['n1','n2'])
	print(tb.shape)
	out_f1 = 'data/ppi/atpin.ppi.2evidence.no-weight.txt'
	out_f2 = 'data/ppi/atpin.ppi.3evidence.no-weight.txt'
	out_f3 = 'data/ppi/atpin.ppi.4evidence.no-weight.txt'
	out_f4 = 'data/ppi/atpin.ppi.1evidence.no-weight.txt'
	out_f5 = 'data/ppi/atpin.ppi.2evidence.weight.txt'
	out_f6 = 'data/ppi/atpin.ppi.3evidence.weight.txt'
	out_f7 = 'data/ppi/atpin.ppi.4evidence.weight.txt'
	out_f8 = 'data/ppi/atpin.ppi.1evidence.weight.txt'
	tb[tb['nu'] > 1][['n1', 'n2', 'w']].to_csv(out_f1, sep=' ', header=None, index=False) 
	tb[tb['nu'] > 2][['n1', 'n2', 'w']].to_csv(out_f2, sep=' ', header=None, index=False) 
	tb[tb['nu'] > 3][['n1', 'n2', 'w']].to_csv(out_f3, sep=' ', header=None, index=False) 
	tb[['n1', 'n2', 'w']].to_csv(out_f4, sep=' ', header=None, index=False) 
	tb[tb['nu'] > 1][['n1', 'n2', 'nu']].to_csv(out_f5, sep=' ', header=None, index=False) 
	tb[tb['nu'] > 2][['n1', 'n2', 'nu']].to_csv(out_f6, sep=' ', header=None, index=False) 
	tb[tb['nu'] > 3][['n1', 'n2', 'nu']].to_csv(out_f7, sep=' ', header=None, index=False) 
	tb[['n1', 'n2', 'nu']].to_csv(out_f8, sep=' ', header=None, index=False) 
def filt_string(f):
	tb = pd.read_csv(f, sep=' ')
	tb.columns = ['n1', 'n2', 'm']
	tb['w'] = 1
	tb = tb.drop_duplicates(subset=['n1','n2'])
	out_f1 = 'data/string.ppi.300.no-weight.txt'		
	out_f2 = 'data/string.ppi.500.no-weight.txt'		
	out_f3 = 'data/string.ppi.700.no-weight.txt'		
	out_f4 = 'data/string.ppi.900.no-weight.txt'		
	out_f5 = 'data/string.ppi.300.weight.txt'		
	out_f6 = 'data/string.ppi.500.weight.txt'		
	out_f7 = 'data/string.ppi.700.weight.txt'		
	out_f8 = 'data/string.ppi.900.weight.txt'
	tb[tb['m'] > 300][['n1', 'n2', 'w']].to_csv(out_f1, sep=' ', header=None, index=False) 
	tb[tb['m'] > 500][['n1', 'n2', 'w']].to_csv(out_f2, sep=' ', header=None, index=False) 
	tb[tb['m'] > 700][['n1', 'n2', 'w']].to_csv(out_f3, sep=' ', header=None, index=False) 
	tb[tb['m'] > 900][['n1', 'n2', 'w']].to_csv(out_f4, sep=' ', header=None, index=False) 
	tb[tb['m'] > 300][['n1', 'n2', 'm']].to_csv(out_f5, sep=' ', header=None, index=False) 
	tb[tb['m'] > 500][['n1', 'n2', 'm']].to_csv(out_f6, sep=' ', header=None, index=False) 
	tb[tb['m'] > 700][['n1', 'n2', 'm']].to_csv(out_f7, sep=' ', header=None, index=False) 
	tb[tb['m'] > 900][['n1', 'n2', 'm']].to_csv(out_f8, sep=' ', header=None, index=False) 
if __name__ == '__main__':
	filt_atpin(ppi_f1)	
#	filt_string(ppi_f2)
	d = 'data/ppi'			
	out_d = 'data/ppi_mult'
	for i in os.listdir(d):
		if i.startswith('atpin.ppi') and i.endswith('txt'):
			mat_f = os.path.join(d, i)
			out = os.path.join(out_d, i)
			sort(mat_f, 'data/arabi_gdic.txt', out=out, rm=None)
		#	sort('data/atpin.ppi.4evidence.no-weight.txt', 'data/gdic.txt', out='data/test_sort.txt', rm=None)
