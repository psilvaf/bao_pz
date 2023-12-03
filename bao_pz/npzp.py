import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
import time
from multiprocessing import Pool


def zp(binsize,z_name,survey_bins,output_folder):
	'''Computes the average number density of the spec-z distribution (not normalised)
	binsize(float): size of photo-z range
	survey_bin(fits Table): data
	output_folder(str): file path
	z_name:name of pz column
	return the number counts of the photometric redshift of a z_dist
	'''	
	COUNT=[np.histogram(survey_bins[l][z_name], bins=binsize) for l in range(len(survey_bins))]
	
	for r in range(len(COUNT)):
		np.save(output_folder+str(r),COUNT[r])
	return

def npdf_zp(bins):
    '''
    Finds the mode of the n(z) values for a bin and computes the representative n(z) from that bin.
    bins(fits/Table): the bins information, where the pdf column must be named 'pdf'
    Return (arr):  the n_p(z_p) for each bin.
    
    '''
    pdf_Count=[[np.histogram(bins[i]['pdf'].T[j]) for j in range(len(bins[i]['pdf'].T))] for i in range(len(bins))]
    maximum=[[] for i in range(len(pdf_Count))]
    for i in range(len(pdf_Count)):
        for j in range(len(pdf_Count[i])):
            maximum[i].append(pdf_Count[i][j][1][np.where(pdf_Count[i][j][0] == max(pdf_Count[i][j][0]))])
    matrix=np.array([np.array([maximum[i][j][0] for i in range(len(maximum))]) for j in range(len(maximum[0]))])
    return np.array(maximum[i]).T[0]/max(np.array(maximum[i]).T[0])


