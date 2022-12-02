import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
import time
from multiprocessing import Pool


def zp(binsize,survey_bins,output_folder):
	'''Computes the average number density of the spec-z distribution (not normalised)
	binsize(float): size of photo-z range
	survey_bin(fits Table): data
	output_folder(str): file path
	
	return the number counts of the photometric redshift of a z_dist
	'''	
	COUNT=[np.histogram(survey_bins[l]['Z'], bins=binsize) for l in range(len(survey_bins))]
	
	for r in range(len(COUNT)):
		np.save(output_folder+str(r),COUNT[r])
	return
