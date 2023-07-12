import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import os
from scipy.signal import find_peaks
from astropy.table import Table
import os.path as path
from scipy import stats

def pdfs_peaks_criteria(data,redshift_name,pdf_name,num_of_peaks):
	'''
	Find pdfs with more peaks that are less than 30% of the main peak
	data: fits table
	redshift_name: redshift columns name
	pdf_name: pdf column name
	num_of_peaks: number of peaks in a pdf
	'''
	peaks=[find_peaks(data[pdf_name][i],height=max(data[pdf_name][i])*.3,)[0] for i in range(len(data))]
	less_peaks=[i for i in range(len(peaks)) if len(peaks[i])<=num_of_peaks]
	return less_peaks
    
def find_multi_modal(data,redshift_name,pdf_name,num_of_peaks):
	'''
	Find pdfs with N peaks 
	data: fits table
	redshift_name: redshift columns name
	pdf_name: pdf column name
	num_of_peaks: number of peaks in a pdf
	'''
	peaks=[find_peaks(data[pdf_name][i])[0] for i in range(len(data))]
	less_peaks=[i for i in range(len(peaks)) if len(peaks[i])==num_of_peaks]	
	return less_peaks
