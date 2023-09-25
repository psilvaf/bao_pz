import numpy as np
from astropy.io import fits
from astropy.table import Table
import os
from os.path import dirname, abspath, join as pjoin
import Corrfunc
from Corrfunc.mocks.DDrppi_mocks import DDrppi_mocks
from astropy.cosmology import FlatLambdaCDM
from Corrfunc.utils import convert_3d_counts_to_cf
from Corrfunc.utils import convert_rp_pi_counts_to_wp
from scipy import interpolate,integrate


def p_corr(survey,random,output_file,weight_name,pimax=120,cosmology=FlatLambdaCDM(H0=67.27, Om0=0.3166, Tcmb0=2.725,Ob0=0.04884)):# planck 18 TT, TE, EE, lowE
	# Setup the bins
	print(cosmology)
	dist=cosmology.comoving_distance(survey['Z'])/0.6727
	#X,Y,Z=dist*survey['mu'],dist*np.sin(np.arccos(survey['mu'])),dist
	print('survey dists')
	dist2=cosmology.comoving_distance(random['Z'])/0.6727
	#X2,Y2,Z2=dist2*random['mu'],dist2*np.sin(np.arccos(random['mu'])),dist2
	print('random dists')
	nthreads = 4
	binfile=np.arange(40,175,5.) #Mpc/h

	autocorr=1

	DD_counts = DDrppi_mocks(autocorr,2, nthreads,pimax, binfile,survey['RA'],survey['DEC'],dist, weights1=survey[weight_name],weight_type='pair_product',output_rpavg=True, is_comoving_dist=True)
	np.save(output_file+'DD',DD_counts)
	RR_counts = DDrppi_mocks(autocorr,2, nthreads,pimax, binfile,random['RA'],random['DEC'],dist2, weights1=random[weight_name],weight_type='pair_product',output_rpavg=True, is_comoving_dist=True)
	np.save(output_file+'RR',RR_counts)
	DR_counts = DDrppi_mocks(0,2, nthreads,pimax, binfile,survey['RA'],survey['DEC'],dist,weights1=survey[weight_name],RA2=random['RA'],DEC2=random['DEC'],CZ2=dist2,weights2=random[weight_name],weight_type='pair_product',is_comoving_dist=True)
	np.save(output_file+'DR',DR_counts)
	return 

def paralell(r,mu): 
	return mu*r 

def perp(r,mu):
	return r*(1-mu**2)**.5

def window(mu):
    if 0<=mu<0.8:
        return 1 #np.random.normal(mu,.3)
    else:
        return 0
        
def interp_xi(perp,par,xi):
	return interpolate.interp2d(perp,par,xi,kind='quintic',fill_value=True)

def function_3d(perp,par,xi_func):
    wp=xi_func(perp,par)
    return wp
    
def spp(mu,perp): 
	return mu*perp/((1-mu**2)**.5)

