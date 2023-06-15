import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
from os.path import dirname, abspath, join as pjoin
import Corrfunc
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from astropy.cosmology import FlatLambdaCDM

def p_corr(survey,random,output_file,weight_name,pimax=175,cosmology=FlatLambdaCDM(H0=70, Om0=0.25, Tcmb0=2.725,Ob0=0.044)):
	# Setup the bins
	print(cosmology)
	dist=cosmology.comoving_distance(survey['Z'])	
	#X,Y,Z=dist*survey['mu'],dist*np.sin(np.arccos(survey['mu'])),dist
	print('survey dists')
	dist2=cosmology.comoving_distance(random['Z'])
	#X2,Y2,Z2=dist2*random['mu'],dist2*np.sin(np.arccos(random['mu'])),dist2
	print('random dists')
	nthreads = 4
	binfile=pjoin(dirname(abspath(Corrfunc.__file__)),"../mocks/tests/", "bins")

	autocorr=1

	DD_counts = DDrppi_mocks(autocorr,2, nthreads,pimax, binfile,survey['RA'],survey['DEC'],dist, weights1=survey[weight_name],weight_type='pair_product',output_rpavg=True, is_comoving_dist=True)
	print('DD')
	RR_counts = DDrppi_mocks(autocorr,2, nthreads,pimax, binfile,random['RA'],random['DEC'],dist2, weights1=random[weight_name],weight_type='pair_product',output_rpavg=True, is_comoving_dist=True)
	print('RR')
	DR_counts = DDrppi_mocks(0,2, nthreads,pimax, binfile,survey['RA'],survey['DEC'],dist,weights1=survey[weight_name],RA2=random['RA'],DEC2=random['DEC'],CZ2=dist2,weights2=random[weight_name],weight_type='pair_product',is_comoving_dist=True)
	
	xi=convert_rp_pi_counts_to_wp(len(survey), len(survey), len(random), len(random),DD_counts['npairs'], DR_counts['npairs'],DR_counts['npairs'], RR_counts['npairs'],14, pimax)
	
	np.save(output_file,xi)
	return xi

