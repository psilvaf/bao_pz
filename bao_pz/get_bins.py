import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import os
import pandas as pd
from multiprocessing import Pool
from astropy.table import Table
import time

def bins(survey,z_name,num,path):
	start = time.time()
	
	data = fits.open(os.path.join(survey))[1].data
	zbin=np.linspace(min(data[z_name]),max(data[z_name]),num)
	with Pool() as pool:
		for i in range(len(zbin)-1):
			tab=data[(zbin[i]<=data[z_name]) & (data[z_name]<zbin[i+1])]
			tabela=Table(tab)
			files=path+str(i)+'.fits'
			tabela.write(files,overwrite=True)
	end = time.time()
	print(end - start)
	return
