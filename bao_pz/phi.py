from astropy.io import fits
import numpy as np
import os
import pandas as pd
from multiprocessing import Pool
from astropy.table import Table
import time

def Phi(z_dist,pdfs_bins,outputfile):    
    z_count=[[[] for i in range(len(z_dist))] for j in range(len(pdfs_bins))]
    with Pool() as pool:
    	for l in range(len(pdfs_bins)):
    		for i in range(len(pdfs_bins[l])):
    			for j in range(len(z_dist)):
    				if not np.interp(z_dist[j],z_dist,pdfs_bins[l][i])==0:
    					z_count[l][j].append(np.interp(z_dist[j],z_dist,pdfs_bins[l][i]))
    	phi0=[[] for i in range(len(z_count))]
    	for i in range(len(z_count)):
    		for j in range(len(z_count[i])):
    			phi0[i].append(len(z_count[i][j]))
    	phi=[i/np.sum(i) for i in phi0]
    	for i in range(len(phi)):
    		np.save(outputfile+str(i),phi[i])
    	return
