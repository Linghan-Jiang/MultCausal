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
warnings.filterwarnings('ignore')
# x = z*bzx + c + ezx
# y = x*bxy + c + exy
def simu_z_bzx(self):
	f = uniform(0.01, 0.5)
	z = binomial(2, f, 100)
	b_zx = normal(0, np.sqrt(1/(2*f*(1-f))))
	z_bzx = np.mean(z*b_zx)
	return z_bzx
def simu_bzx(args):
	mean_z_bzx = args[0]
	var_z_bzx = args[1]
	var_c = args[2]
	var_ezx = args[3]
	f = uniform(0.01, 0.5)
	z = binomial(2, f, 50000).reshape(-1, 1)
	x = normal(mean_z_bzx, np.sqrt(var_z_bzx), 50000) + normal(0, np.sqrt(var_c), 50000) + normal(0, np.sqrt(var_ezx), 50000) 
	x = x.reshape(-1, 1)
	reg = LinearRegression().fit(z, x)
#	reg = LinearRegression(normalize=True).fit(z, x)
	return reg.coef_[0]
def simu(r_zx_2):
	pool = Pool(40)
	r_cx_2 = 0.5	 
	z_bzx = pool.map(simu_z_bzx, range(1000))
	var_z_bzx = np.nanvar(z_bzx)
	mean_z_bzx = np.mean(z_bzx)
	var_c = var_z_bzx*r_cx_2/r_zx_2
	var_ezx = var_z_bzx*(1/(r_zx_2)-1) - var_c
	args = [mean_z_bzx, var_z_bzx, var_c, var_ezx]
	coef = pool.map(simu_bzx, [args for i in range(1000)])
	print(np.mean(coef))	
	print(np.var(coef))

if __name__ == '__main__':
	for r_zx_2 in [0.3, 0.2, 0.15]:
		print('r_zx_2', r_zx_2)
		s = simu(r_zx_2)
