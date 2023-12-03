from astropy.io import fits
import numpy as np
import os
import pandas as pd
from astropy.table import Table
import time
from scipy.interpolate import CubicSpline

def Phi(nz_survey,zp_survey,pdfs_bins,outputfile):
	'''
	Computes the distribution of spec-z given the photo-z of a galaxy sample bin
	with photo-z pdfs
	z_dist(array): redshift distribution expected
	pdfs_bins:list of arrays of the distribution of each object
	output_folder1(str): file path of the counts normalized
	output_folder2(str): file path of the counts averaged
	output_folder3(str): file path of the counts weitght averaged
	
	return: give the distributions normalized 
	'''
	start = time.time()


	z_survey=nz_survey[0]
	z_count=[[[] for i in range(len(z_survey))] for j in range(len(pdfs_bins))]
	for l in range(len(pdfs_bins)):
		for i in range(len(pdfs_bins[l]['pdf'])):
			for j in range(len(z_survey)):
				if not np.interp(z_survey[j],zp_survey,pdfs_bins[l]['pdf'][i])==0:
					z_count[l][j].append(np.interp(z_survey[j],zp_survey,pdfs_bins[l]['pdf'][i]))
	f_w=[[] for i in range(len(z_count))]
	
	for i in range(len(z_count)):
		for j in range(len(z_count[i])):
			f_w[i].append(np.sum(z_count[i][j])*len(z_count[i][j]))
	f_weight=[f_w[i]/np.sum(f_w[i]) for i in range(len(f_w))]
	f_weighted=[f_weight[i]/np.sum(f_weight) for i in range(len(f_weight))]
	for i in range(len(f_weighted)):
		np.save(outputfile+str(i),f_weighted[i])
	end = time.time()
	print(end - start)
	return
	
def wPhi(nz_survey,zp_survey,pdfs_bins,outputfile):
	'''
	Computes the distribution of spec-z given the photo-z of a galaxy sample bin
	with photo-z pdfs
	z_dist(array): redshift distribution expected
	pdfs_bins:list of arrays of the distribution of each object
	output_folder1(str): file path of the counts normalized
	output_folder2(str): file path of the counts averaged
	output_folder3(str): file path of the counts weitght averaged
	
	return: give the distributions normalized 
	'''
	start = time.time()


	z_survey=nz_survey[0]
	
	z_count=[[[] for i in range(len(z_survey))] for j in range(len(pdfs_bins))]
	for l in range(len(pdfs_bins)):
		print('Bin:'+str(l),time.localtime())
		for i in range(len(pdfs_bins[l]['pdf'])):
			ratio=(pdfs_bins[l]['pdf'][i]/pdfs_bins[l]['wfkp_ross'][i])*np.sum(pdfs_bins[l]['wfkp_ross'])
			for j in range(len(z_survey)):
				interpolation=np.interp(z_survey[j],zp_survey,ratio)
				if not interpolation==0:
					z_count[l][j].append(interpolation)
					
					
	f_w=[[] for i in range(len(z_count))]	
	for i in range(len(z_count)):
		for j in range(len(z_count[i])):
			f_w[i].append(np.sum(z_count[i][j])*len(z_count[i][j]))
	f_weight=[f_w[i]/np.sum(f_w[i]) for i in range(len(f_w))]
	f_weighted=[f_weight[i]/np.sum(f_weight) for i in range(len(f_weight))]
	for i in range(len(f_weighted)):
		np.save(outputfile+str(i),f_weighted[i])
	end = time.time()
	print(end - start)
	return

def 
