import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import os
import pandas as pd
from multiprocessing import Pool
from astropy.table import Table
import time
from scipy.interpolate import CubicSpline
from scipy import integrate

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return int(n * multiplier) / multiplier

def Distz(z_dist,zname,binsize,survey,output_folder):
	
	'''Computes the photo-z distribution(matched to the training set) ,
	given the spec-z training set
	z_dist(array): redshift distribution expected
	z_name: name of the redshift column
	binsize(int/array): size of photo-z PDF/defined bins
	survey_bin(fits Table): data
	output_folder(str): file path
	
	return: the number counts of the redshift of a z_dist
	'''
	
	COUNT=[i for i in survey[zname] for j in z_dist if truncate(i,3)==j]
				
				
	np.save(output_folder,np.histogram(COUNT,bins=binsize))
	return 


def NZ_survey(n_train_zs,g_zp,n_train_zp,n_survey):
	'''
	Computes the resampled redshift distribution of the galaxies both in
	the survey and training set. Eq. 4 of arXiv:2210.05057
	n_train_zs: distribution of the spec-z
	g_zp: photo-z distribution(matched to the training set) ,
	given the spec-z training set
	n_train_zp: photo-z of the matched training set
	n_survey: the whole survey photo-z
	'''
	gzp=CubicSpline(g_zp[1][:-1],g_zp[0])
	NZp=CubicSpline(n_train_zp[1][:-1],n_train_zp[0])
	NZ_survey_zp=CubicSpline(n_survey[1][:-1],n_survey[0])
	func=lambda zp: n_train_zs[0]*NZ_survey_zp(zp)*gzp(zp)/NZp(zp)
	
	return n_train_zs[1][:-1],integrate.quad_vec(func,0,2)[0]
	
