from astropy.io import fits
import numpy as np
import os
import pandas as pd
from multiprocessing import Pool
from astropy.table import Table
import time

def Phi(z_dist,pdfs_bins,outputfile1,outputfile2,outputfile3):    
    z_count=[[[] for i in range(len(z_dist))] for j in range(len(pdfs_bins))]
    with Pool() as pool:
    	for l in range(len(pdfs_bins)):
    		for i in range(len(pdfs_bins[l])):
    			for j in range(len(z_dist)):
    				if not np.interp(z_dist[j],z_dist,pdfs_bins[l][i].T[0])==0:
    					z_count[l][j].append(np.interp(z_dist[j],z_dist,pdfs_bins[l][i].T[0]))
    	phi0=[[] for i in range(len(z_count))]
    	f=[[] for i in range(len(z_count))]
    	for i in range(len(z_count)):
    		for j in range(len(z_count[i])):
    			phi0[i].append(len(z_count[i][j]))
    			f[i].append(np.mean(z_count[i][j]))
    			f_w[i].append(np.sum(z_count[i][j])*len(z_count[i][j]))
    	phi=[i/np.sum(i) for i in phi0]
    	f_weighted=[f_w[i]/np.sum(phi0[i]) for i in range(len(phi0))]
    	for i in range(len(phi)):
    		np.save(outputfile1+str(i),phi[i])
    		np.save(outputfile2+str(i),f[i])
    		np.save(outputfile3+str(i),f_weighted[i])
    	return
