import numpy as np
from astropy.io import fits
from astropy.table import Table
import os

def distz(pdfs,output_folder):
	'''Finds the redshfit distribution for each bin
	pdfs: pdf of each galaxy in each bin
	output_folder: name of the output folder, the file will be .npy
	'''
	
	new_dist=[[] for i in range(len(pdfs))]
	for i in range(len(pdfs)):
		for j in range(len(pdfs[i])):
			new_dist[i].append(np.mean(pdfs[i][j],axis=0))
			
	for k in new_dist:
		np.save(output_folder+str(k), k)
	return
	
	
