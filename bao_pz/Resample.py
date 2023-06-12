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
	
def n_z(bins,zname,nz_survey,zp_survey):
	'''
	Computes the average distribution of redshift of each bin, used for both spec-z and photo-z
	bins(dic/astropy table): the separation of the sample in bins
	zname(str): name of the spec-z/photo-z column
	nz_survey(2d array): output from NZ_survey 
	return:  bin name and n_z
	'''
	
	n_z_resampled=CubicSpline(nz_survey[0],nz_survey[1]/np.sum(nz_survey[1]))
	z_survey=nz_survey[0][np.where(n_z_resampled(zp_survey)>10e-3)[0]]
	count=[]
	bin_name=[]
	for i in range(len(bins)):
		count.append([np.where(abs(bins[i][zname]-k)<=.01) for k in z_survey])
		bin_name.append(i)
	normalized_count=[np.array([len(i[0]) for i in j])/np.sum([len(i[0]) for i in j]) for j in count]
	funcs=[CubicSpline(z_survey,np.nan_to_num(normalized_count[i]))(zp_survey) for i in range(len(normalized_count))]
	return bin_name,funcs
	
	


