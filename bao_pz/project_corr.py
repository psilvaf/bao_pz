import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
from os.path import dirname, abspath, join as pjoin
import Corrfunc
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from astropy.cosmology import FlatLambdaCDM

def p_corr(survey,random,output_file,weight_name,pimax=120/.6727,cosmology=FlatLambdaCDM(H0=67.27, Om0=0.3166, Tcmb0=2.725,Ob0=0.04884)):# planck 18 TT, TE, EE, lowE
	# Setup the bins
	print(cosmology)
	dist=cosmology.comoving_distance(survey['Z'])/0.6727
	#X,Y,Z=dist*survey['mu'],dist*np.sin(np.arccos(survey['mu'])),dist
	print('survey dists')
	dist2=cosmology.comoving_distance(random['Z'])/0.6727
	#X2,Y2,Z2=dist2*random['mu'],dist2*np.sin(np.arccos(random['mu'])),dist2
	print('random dists')
	nthreads = 4
	binfile=np.arange(20,175,5.)/.6727 #Mpc/h

	autocorr=1

	DD_counts = DDrppi_mocks(autocorr,2, nthreads,pimax, binfile,survey['RA'],survey['DEC'],dist, weights1=survey[weight_name],weight_type='pair_product',output_rpavg=True, is_comoving_dist=True)
	np.save(output_file+'DD',DD_counts)
	RR_counts = DDrppi_mocks(autocorr,2, nthreads,pimax, binfile,random['RA'],random['DEC'],dist2, weights1=random[weight_name],weight_type='pair_product',output_rpavg=True, is_comoving_dist=True)
	np.save(output_file+'RR',RR_counts)
	DR_counts = DDrppi_mocks(0,2, nthreads,pimax, binfile,survey['RA'],survey['DEC'],dist,weights1=survey[weight_name],RA2=random['RA'],DEC2=random['DEC'],CZ2=dist2,weights2=random[weight_name],weight_type='pair_product',is_comoving_dist=True)
	np.save(output_file+'DR',DR_counts)
	xi=convert_rp_pi_counts_to_wp(len(survey), len(survey), len(random), len(random),DD_counts['npairs'], DR_counts['npairs'],DR_counts['npairs'], RR_counts['npairs'],len(binfile)-1, pimax)
	
	np.save(output_file,xi)
	return 

