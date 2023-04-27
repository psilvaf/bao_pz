from astropy.io import fits
import numpy as np
import os
import pandas as pd
from multiprocessing import Pool
from astropy.table import Table
import time

def Phi(z_dist,survey_bin,outputfile):
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
	
	z_count=[[[] for i in range(len(z_dist))] for j in range(len(survey_bin))]
	with Pool() as pool:
		for l in range(len(survey_bin)):
			for i in range(len(survey_bin[l]['pdf'])):
				for j in range(len(z_dist)):
					if not np.interp(z_dist[j],z_dist,survey_bin[l]['pdf'][i])==0:
						z_count[l][j].append(np.interp(z_dist[j],z_dist,survey_bin[l]['pdf'][i]))
		#counts=[[] for i in range(len(z_count))]
		f=[[] for i in range(len(z_count))]
		f_w=[[] for i in range(len(z_count))]
		for i in range(len(z_count)):
			for j in range(len(z_count[i])):
				#counts[i].append(len(z_count[i][j]))
				f[i].append(np.mean(z_count[i][j]))
				f_w[i].append(np.sum(z_count[i][j])*len(z_count[i][j]))
		#phi=[i/np.sum(i) for i in counts]
		f_weight=[f_w[i]/np.sum(f_w[i]) for i in range(len(f_w))]
		f_weighted=[f_weight[i]/np.sum(f_weight) for i in range(len(f_weight))]
		for i in range(len(f_weighted)):
			#np.save(outputfile1+str(i),phi[i])
			#np.save(outputfile2+str(i),f[i])
			np.save(outputfile+str(i),f_weighted[i])
	end = time.time()
	print(end - start)
	return
