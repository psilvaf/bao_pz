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
	
def moments(pdf,x,n):
    '''pdf(arr): PDf
    n (int): statistical moment
    x(array): variable range
    '''
	mu_n=pdf*(x)**n
	return np.sum(mu_n)
	
def find_gaussians(data,x,output,n=2):

    '''
    Selects the PDFs that are nearly Gaussian.
    
    data (fits loaded file/ astropy Table)
    x
    output (str): output file directory
    n (int): statistical moment
    x(array): variable range
    '''
	pdf=data['pdf']
	mu2_s2=(data['Z_err']**2)+data['Z']**2
	res=np.array([moments(pdf[i],np.linspace(0,2,len(data['pdf'][0])),2) for i in range(len(data))])

	gaussian=np.where((.9<res/mu2_s2) & (res/mu2_s2<=1))[0]
	Table(data[gaussian]).write(output,format='fits',overwrite=True)
	return 


