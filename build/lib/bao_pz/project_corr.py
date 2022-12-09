import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from astropy.cosmology import FlatLambdaCDM

def p_corr(survey,random,output_file,nbins=40, pimax=175,cosmology=FlatLambdaCDM(H0=70, Om0=0.25, Tcmb0=2.725,Ob0=0.044)):
	# Setup the bins
	print(cosmology)
	dist=cosmology.comoving_distance(survey['Z'])
	dist2=cosmology.comoving_distance(random['Z'])
	nthreads = 4
	bins = np.linspace(40, 140.0, nbins + 1)
	pimax=pimax
	autocorr=1
	cosmology=1
	DD_counts = DDrppi_mocks(autocorr, cosmology, nthreads,pimax, bins,survey['RA'],survey['DEC'],dist, survey['WTOT'],weight_type='pair_product',is_comoving_dist=True)
	
	RR_counts = DDrppi_mocks(autocorr, cosmology,  nthreads,pimax, bins,random['RA'], random['DEC'],dist2, random['WEIGHT'],weight_type='pair_product',is_comoving_dist=True)
	
	DR_counts = DDrppi_mocks(autocorr, cosmology, nthreads,pimax, bins,survey['RA'], survey['DEC'],dist,survey['WTOT'],random['RA'], random['DEC'],dist2, random['WEIGHT'],weight_type='pair_product',is_comoving_dist=True)
	
	xi=convert_rp_pi_counts_to_wp(len(survey), len(survey), len(random), len(random),DD_counts['npairs'], DR_counts['npairs'],DR_counts['npairs'], RR_counts['npairs'],nbins, pimax)
	return  np.save(output_file,xi)

