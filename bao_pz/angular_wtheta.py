import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
from Corrfunc.utils import convert_3d_counts_to_cf
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks

def ang_corr(bins,bins2,output_file):
	'''Computes angular correlation function for each bin
	bins: data bins in fits file
	bins2:random bins
	output_file=name of the output file	
	'''
	DDs=[]
	RRs=[]
	DRs=[]
	for i in range(len(bins)):
		# Setup the bins
		nbins = len(bins)
		theta = np.arange(0.5, 5, 22+1) # note the +1 to nbins

		nthreads = 4
		autocorr=1
		DD_counts = DDtheta_mocks(autocorr, nthreads, theta, bins[i]['RA'], bins[i]['DEC'])
		DDs.append(DD_counts)
		
		RR_counts = DDtheta_mocks(autocorr, nthreads, theta, bins2[i]['RA'], bins2[i]['DEC'])
		RRs.append(RR_counts)
		DR_counts = DDtheta_mocks(autocorr, nthreads, theta, bins[i]['RA'], bins[i]['DEC'], RA2=bins2[i]['RA'], DEC2=bins2[i]['DEC'],weights1=bins[i]['Weight'],weights2=bins2[i]['Weight'])
		DRs.append(DR_counts)
		
	wtheta=[]
	for i in range(len(bins)):
		wtheta.append(convert_3d_counts_to_cf(len(bins[i]), len(bins[i]), len(bins2[i]), len(bins2[i]),DDs[i]['npairs'], DRs[i]['npairs'],DRs[i]['npairs'], RRs[i]['npairs']))
	return  np.save(output_file,wtheta)


