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
from colossus.cosmology import cosmology
from multiprocessing import Pool
cosmo=cosmology.setCosmology('planck18')

def opt_dist(x):
    return cosmo.comovingDistance(z_min=x, z_max=0.0, transverse=True)


def DD_mock(survey,output_file,n_threads,pimax=120):
        # Setup the bins

        dist=np.empty(len(survey['Z']))
        with Pool() as p:
                dist=p.map(opt_dist, survey['Z'])

        nthreads = n_threads
        binfile=np.arange(20,175,5.) #Mpc/h

        autocorr=1
        print('Counting DD')
        DD_counts = DDrppi_mocks(autocorr,2, nthreads,pimax,binfile,np.float32(survey['RA']+10),np.float32(survey['DEC']),np.float32(dist), weights1=np.float32([1.]*len(survey)),weight_type='pair_product',output_rpavg=True, is_comoving_dist=True)
        np.save(output_file+'DD',DD_counts)

def p_corr(survey,random,output_file,weight_name,zname,n_threads,pimax=120):# planck 18 TT, TE, EE, lowE
	# Setup the bins
		
	dist=np.empty(len(survey[zname]))
	with Pool() as p:
        	dist=p.map(opt_dist, survey[zname])

	nthreads = n_threads
	binfile=np.arange(20,175,5.) #Mpc/h

	autocorr=1
	print('Counting DD')
	DD_counts = DDrppi_mocks(autocorr,2, nthreads,pimax,binfile,survey['RA'],survey['DEC'],dist, weights1=survey[weight_name],weight_type='pair_product',output_rpavg=True, is_comoving_dist=True)
	np.save(output_file+'DD',DD_counts)
	dist2=np.empty(len(random['Z']))
	with Pool() as p:
		dist2=p.map(opt_dist, random['Z'])

	print('Counting RR')
	RR_counts = DDrppi_mocks(autocorr,2, nthreads,pimax, binfile,random['RA'],random['DEC'],dist2, weights1=random[weight_name],weight_type='pair_product',output_rpavg=True, is_comoving_dist=True)
	np.save(output_file+'RR',RR_counts)
	print('Counting DR')
	DR_counts = DDrppi_mocks(0,2, nthreads,pimax, binfile,survey['RA'],survey['DEC'],dist,weights1=survey[weight_name],RA2=random['RA'],DEC2=random['DEC'],CZ2=dist2,weights2=random[weight_name],weight_type='pair_product',is_comoving_dist=True)
	np.save(output_file+'DR',DR_counts)
	return 

def paralell(r,mu): 
	return mu*r 

def perp(r,mu):
	return r*(1-mu**2)**.5

def window(x):
    if 0<x<20/.6:
        return 1
    else:
        return 0.0

def function_3d(perp,par,xi_func):
    wp=xi_func(perp,par)
    return wp
    
def spp(mu,perp): 
	return mu*perp/((1-mu**2)**.5)
	
def get_func(ND,NR,DD,DR,RR,sperp, spar,n):
    '''Computes the projected correlation function numerically by integrating over spar
    ND(int) : number of galaxies
    NR (int): number of random
    DD(arr): pair count
    DR(arr): pair count
    RR(arr): pair count
    sperp(arr): distance perp to the LOS
    spar(arr): distance par to the LOS
    
    return: projected 2pcf normalized.
    '''
    dd = np.empty(len(sperp))
    dr=np.empty(len(sperp))
    rr=np.empty(len(sperp))
    for i in range(len(sperp)):
        dd[i] = np.sum(DD[i * len(spar):(i + 1) * len(spar)])/ND
        dr[i] = np.sum(DR[i * len(spar):(i + 1) * len(spar)])/(ND*NR)
        rr[i] = np.sum(RR[i * len(spar):(i + 1) * len(spar)])/NR
    corr=(dd-2*dr+rr)/(rr)
    sigma=np.mean(np.round(corr,n))
    return corr-sigma
    
def int_from_mu(perp,r,mu,xi,perp_size,par_size):
	'''Computes the projected correlation function numerically by integrating over mu
	
	perp(arr): chosen perp separation 
	r(arr): full separation in 3D
	mu(arr): mu
	xi(function): correlation function in 3D
	perp_size(int): size of perp separation
	par_size(int): size of par separation
	
	return: projected 2pcf (not normalised)	
	'''
	g= np.array([window(mu[i])*function_3d(perp,paralell(r,mu[i]),xi) for i in range(len(mu))]).sum(axis=0).sum(axis=1)
	G=[np.sum(g[i * par_size:(i + 1) * par_size]) for i in range(perp_size)]
	return G


