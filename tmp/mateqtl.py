import os
from multiprocessing import Pool
cis = 1e-2
tra = 1e-2
pool = Pool(10)
for dis in range(5, 8):
	dis = 10**dis
	pool.apply_async(os.system, ('Rscript mateqtl.R {} {} {}'.format(cis, tra, dis),))
pool.close()
pool.join()		
