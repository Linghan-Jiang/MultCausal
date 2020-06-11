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
def z_bzx_1():
	# pleiotropy model
	# x = zc*bzcx + zp * bzpx + c + ezx
	# y = x*bxy + zp*bzpy + c + exy
	f = uniform(0.01, 0.5)
	z = binomial(2, f, 100)
	b_zx = normal(0, 1/(2*f*(1-f)))
	z_bzx = sum(z*b_zx)
	return z_bzx
def bzx_1(args):
	mean_z_bzx = args[0]
	var_z_bzx = args[1]
	var_c = args[2]
	var_ezx = args[3]
	f = uniform(0.01, 0.5)
	z = binomial(2, f, 50000).reshape(-1, 1)
	x = []
	for i in range(50000):
		x.append(normal(mean_z_bzx, var_z_bzx) + normal(0, var_c) + normal(0, var_ezx)) 
	x = np.array(x).reshape(-1, 1)
	#print(z.shape)
	#print(x.shape) 
	reg = LinearRegression().fit(z, x)
	return reg.coef_
def bxy_1(args):
	mean_z_bzx = args[0]
	var_z_bzx = args[1]
	var_c = args[2]
	var_ezx = args[3]
	var_ey = args[4]
	x = []
	y = []		
	for i in range(20000):
		y.append(normal(0, var_c) + normal(0, var_ey))
		x.append(normal(mean_z_bzx, var_z_bzx) + normal(0, var_c) + normal(0, var_ezx)) 
	y = np.array(y).reshape(-1, 1)
	x = np.array(x).reshape(-1, 1)
	reg = LinearRegression().fit(y, x)			
	return reg.coef_	
def bxy_2(args):
	mean_z_bzx = args[0]
	var_z_bzx = args[1]
	var_c = args[2]
	var_ezx = args[3]
	mean_x_bxy = args[4]
	var_x_bxy = args[5]
	var_exy = args[6]
	
	x = []
	y = []		
	for i in range(20000):
		y.append(normal(mean_x_bxy, var_x_bxy) + normal(0, var_c) + normal(0, var_exy))	
		x.append(normal(mean_z_bzx, var_z_bzx) + normal(0, var_c) + normal(0, var_ezx)) 
	y = np.array(y).reshape(-1, 1)
	x = np.array(x).reshape(-1, 1)
	reg = LinearRegression().fit(y, x)			
	return reg.coef_	
	
class simu_1(object):
	def __init__(self, r_zx_2):
		r_cx_2 = 0.5	 
		z_bzx = []
		for i in range(1000):
			z_bzx.append(z_bzx_1())
		var_z_bzx = np.nanvar(z_bzx)
		mean_z_bzx = np.mean(z_bzx)
		var_c = var_z_bzx*r_cx_2/r_zx_2
		var_ezx = var_z_bzx*(1/(r_zx_2)-1) - var_c
		args = [mean_z_bzx, var_z_bzx, var_c, var_ezx]
		self.mean_z_bzx = mean_z_bzx
		self.var_z_bzx = var_z_bzx
		self.var_ezx = var_ezx
		self.var_c = var_c
		self.r_cx_2 = r_cx_2
	#	coef = pool.map(bzx_1, [args for i in range(1000)])
	#	print(np.mean(coef))
	def worker(self):
		pool = Pool(30)
		r_cy_2 = 0.5
		var_ey = self.var_c*(1/r_cy_2-1)
		self.args.append(var_ey)
		coef = pool.map(bxy_1, [self.args for i in range(1000)])
		print(np.mean(coef))
		print(np.nanvar(coef))	
def x_bxy(args):
	mean_z_bzx = args[0]
	var_z_bzx = args[1]
	var_c = args[2]
	var_ezx = args[3]
	x = []
	for i in range(1000):
		x.append(normal(mean_z_bzx, var_z_bzx) + normal(0, var_c) + normal(    0, var_ezx))
	var_x_bxy = np.nanvar(x)
	mean_x_bxy = np.mean(x)
	return (var_x_bxy, mean_x_bxy)
class simu_2(simu_1):
	def __init__(self, r_zx_2):
		simu_1.__init__(self, r_zx_2)
	def worker(self):
		pool = Pool(40)
		for r_xy_2 in [0.05, 0.1]:
			#r_xy_2 = 0.05
			var_c = self.var_z_bzx * self.r_cx_2/r_zx_2
			args = [self.mean_z_bzx, self.var_z_bzx, self.var_c, self.var_ezx]
			for i in range(1000):
				var_x_bxy, mean_x_bxy = x_bxy(args)
			var_exy = var_x_bxy*(1/r_xy_2-1) - var_c*3
			args = args + [mean_x_bxy, var_x_bxy, var_exy]
			print(args)
			coef = pool.map(bxy_2, [args for i in range(1000)])
			print(np.mean(coef))
			print(np.nanvar(coef))

if __name__ == '__main__':
	for r_zx_2 in [0.15, 0.2, 0.3]:
		#s = simu_1(r_zx_2)
		#s.worker()
		s2 = simu_2(r_zx_2)
		s2.worker()
#	simu_1_worker()
#	simu_1_worker(0)
	
