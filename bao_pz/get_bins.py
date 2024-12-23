import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import os
import pandas as pd
from multiprocessing import Pool
from astropy.table import Table
import time

def bins(survey,z_name,path,zp1,zp2,width):
    '''
    Divides the survey in redhsift bins
    survey (str): fits file name
    z_name (str): name of the photo-z  column name
    path (str): name of the output fits files
    zp1 (float): lowest redshift
    zp2 (float): highest redshift
    width (float): bin separation
    
    return:
    fits files for each bin
    '''



	start = time.time()
	
	data = fits.open(os.path.join(survey))[1].data
	zbin=np.arange(zp1,zp2,width)
	print(zbin)
	with Pool() as pool:
		for i in range(len(zbin)-1):
			tab=data[(data[z_name]>zbin[i]) & (data[z_name]<=zbin[i+1])]
			tabela=Table(tab)
			files=path+str(i)+'.fits'
			tabela.write(files,overwrite=True)
	end = time.time()
	print(end - start)
	return
