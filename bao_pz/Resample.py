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
	
	COUNT=np.histogram(survey[zname][np.where(abs(z_dist-survey[zname])<=.01)],bins=binsize)
				
				
	np.save(output_folder,COUNT)
	return 
	
def NZ_survey(n_train_zs,g_zp,n_train_zp,n_survey):
	'''
	Computes the resampled spec-z distribution of the surbey
	n_train_zs(2d-array): distribution of the training set by a chosen number of separation bins
	g_zp(2d-array): output of Distz
	n_train_zp(2d-array): the photo-z of the training set
	n_survey(2d-array): histogram output of the whole photo-z survey with the same bin separation as n_train_zs and n_train_zp
	return: spec_z_survey,NZ_survey
	'''
	gzp=CubicSpline(g_zp[1][:-1],g_zp[0]/np.sum(g_zp[0]))
	NZp=CubicSpline(n_train_zp[1][:-1],n_train_zp[0]/np.sum(n_train_zp[0]))
	NZ_survey_zp=CubicSpline(n_survey[1][:-1],n_survey[0]/np.sum(n_survey[0]))
	func=lambda zp: n_train_zs[0]*NZ_survey_zp(zp)*gzp(zp)/NZp(zp)
	
	return n_train_zs[1][:-1],integrate.quad_vec(func,0,np.inf)[0]
	
def n_z(bins,zname,nz_survey):
	count=[]
	for i in range(len(bins)):
		count.append([np.where(abs(bins[i][zname]-k)<=.001) for k in nz_survey])

	normalized_count=[np.array([len(i[0]) for i in j])/np.sum([len(i[0]) for i in j]) for j in count]
	return normalized_count
	


