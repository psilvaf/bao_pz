import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
import time
from multiprocessing import Pool

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

def distz(z_dist,binsize,survey_bins,output_folder):
	
	'''Computes the average number density of the spec-z distribution (not normalised)
	z_dist(array): redshift distribution expected
	binsize(int/array): size of photo-z PDF/defined bins
	survey_bin(fits Table): data
	output_folder(str): file path
	
	return: the number counts of the redshift of a z_dist
	'''
	
	COUNT=[[] for j in range(len(survey_bins))]
	for l in range(len(survey_bins)):
		for i in survey_bins[l]['Z']:
			for j in z_dist:
				if truncate(i,2)==j:
					COUNT[l].append(truncate(i,2))					
				
	for r in range(len(COUNT)):
		np.save(output_folder+str(r),np.histogram(COUNT[r],bins=binsize))
	return


def distzp(binsize,survey_bins,output_folder):
	'''Computes the average number density of the spec-z distribution (not normalised)
	binsize(float): size of photo-z PDF/defined bins
	survey_bin(fits Table): data
	output_folder(str): file path
	
<<<<<<< Updated upstream
	return the number counts of the photometric redshift of a z_dist
	'''	
	COUNT=[np.histogram(survey_bins[l]['Z'], bins=binsize) for l in range(len(survey_bins))]
	
=======
def distzp(binsize,survey_bins,output_folder):
	'''Computes the average number density of the spec-z distribution (not normalised)
	binsize(float): size of photo-z PDF/defined bins
	survey_bin(fits Table): data
	output_folder(str): file path
	
	return the number counts of the photometric redshift of a z_dist
	'''	
	COUNT=[np.histogram(survey_bins[l]['Z'], bins=binsize) for l in range(len(survey_bins))]
	
>>>>>>> Stashed changes
	for r in range(len(COUNT)):
		np.save(output_folder+str(r),COUNT[r])
	return
