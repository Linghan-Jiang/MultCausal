import os
import numpy as np
from numpy.random import *
import pandas as pd
import networkx as nx
import warnings
#from scipy.stats import 
from matplotlib import pyplot as plt
from multiprocessing import Pool
from utils import num_2_ped
from sklearn.linear_model import LinearRegression
from sklearn import preprocessing
warnings.filterwarnings('ignore')
# x = z*bzx + c + ezx
# y = x*bxy + c + exy
def simu_z_bzx(self=True):
	f = uniform(0.01, 0.5)
	z = binomial(2, f, 100)
	b_zx = normal(0, np.sqrt(1/(2*f*(1-f))))
	z_bzx = np.mean(z*b_zx)
	return z_bzx
def simu_bxy(args):
	place = 0
	mean_z_bzx = args[0]
	var_z_bzx = args[1]
	var_c = args[2]
	var_ezx = args[3]
	mean_x_bxy = args[4]
	var_x_bxy = args[5]
	var_exy = args[6]
#	print('mean_x_bxy', mean_x_bxy)
#	print('mean_z_bzx', mean_z_bzx)
#	x = np.array([simu_z_bzx() for i in range(20000)])
	pool = Pool(40)
	x = pool.map(simu_z_bzx, [1]*20000)
	
	x = x + normal(0, np.sqrt(var_c), 20000) + normal(0, np.sqrt(var_ezx), 20000)
	y = x + normal(0, np.sqrt(var_c), 20000) + normal(0, np.sqrt(var_exy), 20000)
#	x = preprocessing.scale(x, axis=1) 
#	y = preprocessing.scale(y, axis=1)	
#	print(x)
#	with open('xy', 'w+') as f:
#		for i,j in zip(list(x), list(y)):
#			f.write(str(i)+' '+str(j)+'\n')
#	raise
	y = y.reshape(-1, 1)
	x = x.reshape(-1, 1)
#	print('var_y', np.var(y))
#	print('mean_y', np.mean(y))
#	reg = LinearRegression(normalize=True).fit(x, y)			
	reg = LinearRegression().fit(x, y)			
#	print(reg.coef_[0])
	return reg.coef_[0]	
def simu(r_zx_2):
	pool = Pool(2)
	# z_bzx
	r_cx_2 = 0.5	 
	z_bzx = pool.map(simu_z_bzx, range(1000))
	var_z_bzx = np.nanvar(z_bzx)
	mean_z_bzx = np.mean(z_bzx)
	var_c = var_z_bzx*r_cx_2/r_zx_2
	var_ezx = var_z_bzx*(1/(r_zx_2)-1) - var_c
	x = normal(mean_z_bzx, np.sqrt(var_z_bzx), 1000) + normal(0, np.sqrt(var_c), 1000) + normal(0, np.sqrt(var_ezx), 1000)
	var_x_bxy = np.nanvar(x)
	mean_x_bxy = np.mean(x)
	print('var_x_bxy', var_x_bxy)
	print('var_c', var_c)
	print('var_ezx', var_ezx)
	print('var_z_bzx', var_z_bzx)
	print('mean_x_bxy', mean_x_bxy)
	print('mean_z_bzx', mean_z_bzx)
	args = [mean_z_bzx, var_z_bzx, var_c, var_ezx, mean_x_bxy, var_x_bxy]
	for r_xy_2 in [0.05, 0.1]:
		var_exy = var_x_bxy*(1/r_xy_2-1) - var_c*3
		tmp_args = args + [var_exy]
#		simu_bxy(tmp_args)
		coef = pool.map(simu_bxy, [tmp_args for i in range(1000)])
		print('mean bxy', np.mean(coef))
		print('std_bxy', np.std(coef))

if __name__ == '__main__':
	for r_zx_2 in [0.3, 0.2, 0.15]:
		s = simu(r_zx_2)
