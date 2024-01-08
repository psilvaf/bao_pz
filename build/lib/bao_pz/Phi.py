from astropy.io import fits
import numpy as np
import os
import pandas as pd
from astropy.table import Table
import time
from scipy.interpolate import CubicSpline

def compute_phi_numerical(f,nz,npzp,):
	'''
	Computing phi(z|zp)
	f (arr) : the distribution f(z|zp)
	nz (arr) : spec-z distribution of the training set
	npzp (arr) : photo-z distribution of the training set
	
	return (arr) : phi(z|zp)
	'''
	
	n=(nz/np.sum(nz))/(npzp/np.sum(npzp))
	phi=f*n
	return phi
	
def compute_phi(f,nz,npzp,z):
	'''
	Computes phi(z|zp) function
	f (arr) : the distribution f(z|zp)
	nz (arr) : spec-z distribution of the training set
	npzp (arr) : photo-z distribution of the training set
	z (arr) : the spec-z distribution from the resampled distribution
	
	return (function) : phi(z|zp)
	'''
	
	f_func=CubicSpline(z,f)
	nz_func=CubicSpline(z,nz/np.max(nz))
	npzp_func=CubicSpline(z,npzp/np.max(npzp))
	phi=lambda z: f_func(z)*nz_func(z)/npzp_func(z)
	return phi
