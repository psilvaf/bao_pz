import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from astropy.cosmology import FlatLambdaCDM

def ang_corr(survey,random,binfile,cosmology=FlatLambdaCDM(H0=70, Om0=0.25, Tcmb0=2.725,Ode0=0.75,Ob0=0.044),output_file):
	# Setup the bins

	dist=cosmology.comoving_distance(survey['Z'])
	dist2=cosmology.comoving_distance(random['Z'])
	nthreads = 4
	pimax=120
	DD_counts = DDrppi_mocks(autocorr=1, cosmology=1, nthreads,pimax, binfile,RA1=survey['RA'],DEC1= survey['DEC'],CZ1=dist, weights1=survey['WTOT'],weight_type='pair_product',is_comoving_dist=True)

	
	RR_counts = DDrppi_mocks(autocorr=1, cosmology=1, nthreads,pimax, binfile,RA1=random['RA'], DEC1=random['DEC'],CZ1=dist2, weights1=random['WEIGHT'],weight_type='pair_product',is_comoving_dist=True)

	DR_counts = DDrppi_mocks(autocorr=1, cosmology=1, nthreads,pimax, binfile,RA1=survey['RA'], DEC1=survey['DEC'],CZ1=dist,weights1=survey['WTOT'],RA2=random['RA'], DEC2=random['DEC'],CZ2=dist2, weights2=random['WEIGHT'],weight_type='pair_product',is_comoving_dist=True)
		
	xi=convert_rp_pi_counts_to_wp(len(survey), len(survey), len(random), len(random),DD_counts['npairs'], DR_counts['npairs'],DR_counts['npairs'], RR_counts['npairs'])
	
	
	return  np.save(output_file,xi)

