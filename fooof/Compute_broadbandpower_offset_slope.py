#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Computes broadband power, offset and slope of power spectrum
Created on Mon Sep 14 12:09:27 2020

Based on selected epochs (e.g. ASCIIS) in the NO-cohorten file a power spectrum 
is computed. Based on this power spectrum the broadband power is calculated, 
followed by the offset and slope using the FOOOF algorithm. 

reference paper FOOOF: Haller M, Donoghue T, Peterson E, Varma P, Sebastian P, Gao R, Noto T, 
Knight RT, Shestyuk A, Voytek B (2018) Parameterizing Neural Power Spectra. 
bioRxiv, 299859. doi: https://doi.org/10.1101/299859
reference github: https://fooof-tools.github.io/fooof/index.html

@author: t.numan
"""

__author__ = "Tianne Numan"
__contact__ = "t.numan@amsterdamumc.nl" # or l.douw@amsterdamumc.nl
__date__ = "2020/09/14"   ### Date it was created
__status__ = "Production"


####################
# Review History   #
####################

# Reviewed by Name Date ### e.g. Eduarda Centeno 20200909


####################
# Libraries        #
####################

# Standard imports  ### (Put here built-in libraries - https://docs.python.org/3/library/)
import time
import os
import glob

# Third party imports ### (Put here third-party libraries e.g. pandas, numpy)
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import signal
from fooof import FOOOF


# ------------------------------------------------------------------------------
def get_names_epochs(location, nr_epochs=10, extension='.asc'):
    """ obtains the names of MEG ASCIIs as saved in NO-cohorten
    Parameters
    ----------
    location : str
        give the directory where the ASCIIS are stored 
    nr_epochs : int, optional 
        number of epochs, default=10
    extension : str, optional
        give the extension of the file you're looking for, default='.asc'    
    
    Returns
    -------
    selected_asciis : list
        list with strings of ASCII names
    """
    
    selected_asciis = []
    if str(glob.glob(location + '*' + 'Tr_01' + extension)).strip('[]') != "":  # if first case has Tr_01
        for i in range(nr_epochs):
            selected_asciis.append(str(glob.glob(location + '*' + 'Tr_' + "%02d" % (i+1) + extension)).strip('[]'))
    elif str(glob.glob(location + '*' + 'Tr_1' + extension)).strip('[]') != "":  # if first case has Tr_1
        for i in range(nr_epochs):
            selected_asciis.append(str(glob.glob(location + '*' + 'Tr_' + str(i+1) + extension)).strip('[]'))
    elif str(glob.glob(location + '*' + 'Tr01' + extension)).strip('[]') != "":  # if first case has Tr01
        for i in range(nr_epochs):
            selected_asciis.append(str(glob.glob(location + '*' + 'Tr' + "%02d" % (i+1) + extension)).strip('[]'))
    elif str(glob.glob(location + '*' + 'Tr1' + extension)).strip('[]') != "":  # if first case has Tr1
        for i in range(nr_epochs):
            selected_asciis.append(str(glob.glob(location + '*' + 'Tr' + str(i+1) + extension)).strip('[]'))
    return selected_asciis            


def cal_power_spectrum (timeseries, nr_rois=np.arange(92), fs=1250, window='hamming', nperseg=4096, 
                        scaling='spectrum', plot_figure=False, title_plot='average power spectrum'):
    """ calculate (and plot) power spectrum of timeseries
    Parameters
    ----------
    timeseries : DataFrame with ndarrays
        rows are timepoints, columns are rois/electrodes
    rois: int, optional
        give list with rois/electrodes you want to include, default=np.arange(92)
    fs: int, optional    
        sample freuqency, default=1250    
    window: str or tuple, optional
        type of window you want to use, check spectral.py for details, default=' hamming'
    nperseg : int, optional    
        length of each segment, default=4096
    scaling : str, optional
        'density' calculates the power spectral density (V**2/Hz), 'spectrum' calculates the 
        power spectrum (V**2), default='spectrum'
    plot_figure : bool
        creates a figure of the mean + std over all rois/electrodes, default=False
    title_plot : str
        give title of the plot, default='average power spectrum'
        
    Returns
    -------
    f : ndarray
        array with sample frequencies (x-asix of power spectrum plot)
    pxx : ndarray
        columns of power spectra for each roi/VE
            
    """
  
    pxx = np.empty([int(nperseg/2+1), np.size(nr_rois)])
    
    i=0
    for roi in nr_rois:
        (f, pxx[:,i]) = signal.welch(timeseries[roi].values, fs, window, nperseg, scaling=scaling)
        i=i+1
    if plot_figure==True:
        plt.figure()
        plt.plot(f, np.mean(pxx,1), color='teal')
        plt.plot(f, np.mean(pxx,1)+np.std(pxx,1), color='teal', linewidth=0.7)
        plt.plot(f, np.mean(pxx,1)-np.std(pxx,1), color='teal', linewidth=0.7)
        plt.fill_between(f, np.mean(pxx,1)+np.std(pxx,1), np.mean(pxx,1)-np.std(pxx,1), color='teal', alpha=0.2)
        plt.xlim(0, 50)
        plt.xlabel('Frequency (Hz)')
        plt.title(title_plot)
        plt.show()
    return f, pxx


def find_nearest(array, value):
    """ find nearest value of interest in array (used for frequencies, no double value issues)
    Parameters
    ----------
    array : array
        give the array in which you want to find index of value nearest-by
    value : int or float
        the value of interest
     
    Return
    ------
    idx : int
        index of value nearest by value of interest    
    
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def run_loop_powerspectrum(subject_list, extension='.asc', nr_epochs=10, nr_rois=np.arange(92), 
                           Fs=1250, window_length=4096, freq_range=[0.5, 48]):
    """ calculate power spectrum for all cases within the subject_list
    Parameters
    ----------
    subject_list : DataFrame or list
        DataFrame should contain a column 'location_MEG' with the directory of all MEG 
        ASCIIs for each subject. List should contain only 1 column describing the location. 
    extension : str, optional
        give the extension of ASCIIs, '.txt' or '.asc', default='.asc'
    nr_rpochs : int, optional
        give number of epochs you want to include, default=10
    nr_rois : int, optional
        give list with  rois you want to analyse, default=np.arange(92)
    Fs : int, optional
        sample frequency, default=1250
    window_length : int, optional
        window length to calculate power spectrum, default=4096
    freq_range : list, optional
        gives the upper and lower boundaries to calculate the broadband power, default=[0.5, 48]    
    
    Return
    ------
    mean_pxx : ndarray (size: len(subjects), len(power spectrum), nr_rois)
        power spectum for all subjects, and all rois/VE, averaged over nr_epochs
    broadband_power : ndarray (size: len(subjects), nr_rois)
        broadband power between freq_range
    f : ndarray
        array with sample frequencies (x-asix of power spectrum plot)
   
    """  
    
    mean_pxx = np.empty([len(subject_list), int(window_length/2+1), np.size(nr_rois)])
    broadband_power = np.empty([len(subject_list), np.size(nr_rois)])
    if ~isinstance(subject_list, pd.DataFrame): # check if subject_list is dataframe (or list)
        subject_list=pd.DataFrame(subject_list, columns = ['location_MEG']) # save dataframe with column 'location'
    for sub in range(len(subject_list)):
        asciis = get_names_epochs(subject_list['location_MEG'][sub], nr_epochs, extension)
        sub_pxx = np.zeros((nr_epochs, int(window_length/2+1), np.size(nr_rois)))
        mean_pxx[sub,:,:] = 'nan'
        broadband_power[sub,:] = 'nan'
        print(sub)
        if len(asciis) == 0:
            print("no ASCIIs available for: " + subject_list['location'][sub])
            continue    
        for j in range(nr_epochs):
            timeseries =pd.read_csv(asciis[j][1:-1], index_col=False, header=None, delimiter='\t') 
            # first and last character should be removed otherwise it's not working)
            # compute power spectrum
            f, pxx = cal_power_spectrum(timeseries, nr_rois=nr_rois, fs=Fs, plot_figure=False, 
                                        title_plot='average power spectrum - epoch: '+ str(j))
            sub_pxx[j,:,:] = pxx           
        mean_pxx[sub,:,:] = np.nanmean(sub_pxx, axis=0)    
        broadband_power[sub,:] = np.sum(mean_pxx[sub,find_nearest(f,freq_range[0]):find_nearest(f,freq_range[1]),:], 
                                        axis=0)
    return mean_pxx, broadband_power, f  


def cal_FOOOF_parameters(pxx, f, freq_range=[0.5, 48]):
    """ obtain slope and offset using the FOOOF algorithm
    reference paper: Haller M, Donoghue T, Peterson E, Varma P, Sebastian P, Gao R, Noto T, 
    Knight RT, Shestyuk A, Voytek B (2018) Parameterizing Neural Power Spectra. 
    bioRxiv, 299859. doi: https://doi.org/10.1101/299859
    
    reference github: https://fooof-tools.github.io/fooof/index.html
    
    Parameters
    ----------
    pxx : ndarray
        column of power spectra 
    f : 1d-array
        array with sample frequencies (x-asix of power spectrum plot)
    freq_range : list, optional
        gives the upper and lower boundaries to calculate the broadband power, default=[0.5, 48]
    
    Returns
    -------
    FOOOF_offset : float
        offset 
    FOOOF_slope  : float
        slope  
    
    """
    
    # initialize FOOOF oject
    fm = FOOOF()
    # create variables
    fm.fit(f, pxx, freq_range)
    FOOOF_offset = fm.background_params_[0]
    FOOOF_slope = fm.background_params_[1]
    time.sleep(1) # heavy algorithm
    
    return FOOOF_offset, FOOOF_slope


# ------------------------------------------------------------------------------

# set nice level to 10, especially FOOOF algorithm is heavy!
os.nice(10)


###########################
# Settings                #   
###########################

atlas = 'AAL' # AAL or BNA
# select either a list of rois (e.g. np.arange(78) when using the cortical regions of the AAL)
# or a single roi (e.g. (10,) if you want to evaluate roi 11
nr_rois = np.arange(78) # give list of rois or 1 rois (AAL cortical = 78, BNA cortical = 210)
#nr_rois = (10,) # only run for roi 10 # note that python idexes at 0!

nr_epochs= 10 # number of epochs
Fs = 1250 # sample frequency
freq_range=[0.5, 48] # frequency range you want to analyze

subject_list = '/data/KNW/t.numan/GOALS/example_MEG_list.csv'
dir_input = '/data/KNW/NO-cohorten/Scans/'
dir_output = '/data/KNW/t.numan/GOALS/tmp/'
save_output = True # you can save output
name_to_save = 'test' # provide relevant name to save your run


###########################
# Run analysis            # 
###########################

# read the csv file containing Case_ID and MM of each included MEG
# see example '/data/KNW/t.numan/GOALS/example_MEG_list.csv'
subjects = pd.read_csv(subject_list, delimiter=';', header=0, usecols=[0,1])
subjects['Case_ID'] = subjects['Case_ID'].apply(lambda x: str(x).zfill(4)) # add padding zeros when necessary

# get location of ASCIIs of MEG 
subjects['location_MEG'] = dir_input + 'sub-'+ subjects['Case_ID'] + '/meg/' + subjects['MM'] + '/' + atlas + '/' 


###### Compute power spectrum ######
## takes a while when multiple subjects are analyzed at once
mean_pxx, broadband_power, f = run_loop_powerspectrum(subjects, nr_epochs=nr_epochs, nr_rois=nr_rois, 
                                                      Fs=Fs, window_length=4096, freq_range=freq_range)
       
# mean_pxx contains the average power spectra over nr_epochs for all subjects for all rois
# broadband_power is the area under the average power spectra over the frequency range for all subjects and rois
# f gives the frequencies of the power spectrum (can be useful when plotting power spectrum)

# save output
if save_output == True:
    np.save(dir_output + name_to_save + '_mean_pxx', mean_pxx)
    np.save(dir_output + name_to_save + '_broadband_power', broadband_power)
    print('saving power spectra and broadband power')
    



###### Compute FOOOOF offset and slope ######
# create empty arrays to store offset and slope values
offset = np.empty([mean_pxx.shape[0], mean_pxx.shape[2]]) # index 0 -> number of subjects, index 2 -> number of rois
slope = np.empty([mean_pxx.shape[0], mean_pxx.shape[2]])

print('Running FOOOF.........................')

# run across all subjects in your list
for sub in range(len(subjects)):
    print('subject : ' + str(sub)) # print which subject is analyzed
    if np.isnan(mean_pxx[sub,1,0]): # if there is no mean_pxx present for a subject, set offset/slope to nan
        offset[sub,] = np.nan
        slope[sub,] = np.nan
    else:
        i=0
        for roi in nr_rois:
            print('roi : ' + str(roi))
            offset[sub,i], slope[sub,i] = cal_FOOOF_parameters(mean_pxx[sub,:,i], f, freq_range=[0.5, 48])
            i=i+1
            time.sleep(0.05) # by having a pause in the analyses it should limit the %CPU continuously
    # might be relevant when computing FOOOF for large dataset        
    if sub % 10 ==0 & save_output==True:
        print('save offset and slope in between after (another) 10 subjects')
        np.save(dir_output + name_to_save + '_FOOOF_offset', offset)
        np.save(dir_output + name_to_save + '_FOOOF_slope', slope)
    
if save_output == True:          
    np.save(dir_output + name_to_save + '_FOOOF_offset', offset)
    np.save(dir_output + name_to_save + '_FOOOF_slope', slope)
    # save offset and slope output of FOOOF
    print("FOOOF-offset is saved: "+ dir_output + name_to_save + "_FOOOF_offset" )
    print("FOOOF-slope is saved: "+ dir_output + name_to_save + "_FOOOF_slope" )
    
######## Combine results with Case_ID in 'subjects' variable
    
subjects['mean_broadbandpower'] = np.mean(broadband_power, axis=1)     
subjects['mean_offset'] = np.mean(offset, axis=1)
subjects['mean_slope'] = np.mean(slope, axis=1)

subjects.to_csv(dir_output + name_to_save + '_summary.csv' , sep=';', index=False)
# you can import this .csv file in SPSS or Excel

