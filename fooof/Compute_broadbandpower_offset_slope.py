#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Computes broadband power, offset and slope of power spectrum

Based on selected epochs (e.g. ASCIIS) in the NO-cohorten file a power spectrum 
is computed. Based on this power spectrum the broadband power is calculated, 
followed by the offset and slope using the FOOOF algorithm. 

Reference paper FOOOF: Haller M, Donoghue T, Peterson E, Varma P, Sebastian P, 
Gao R, Noto T, Knight RT, Shestyuk A, Voytek B (2018) Parameterizing Neural 
Power Spectra. bioRxiv, 299859. doi: https://doi.org/10.1101/299859
reference Github: https://fooof-tools.github.io/fooof/index.html

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

# Standard imports  
import time
import os
import glob

# Third party imports 
import numpy as np # version 1.19.1
import matplotlib.pyplot as plt # version 3.3.0
import pandas as pd # version 1.1.0
from scipy import signal # version 1.4.1
from fooof import FOOOF # version 0.1.3


# -----------------------------------------------------------------------------
# Define Functions 
def get_names_epochs(location, nr_epochs=10, extension='.asc'):
    """ Obtains the names of MEG ASCIIs as saved in folder (e.g., NO-cohorten)
    
    Parameters
    ----------
    location: str
        Give the directory where the ASCIIS are stored 
    nr_epochs: int, optional 
        Number of epochs, default=10
    extension: str, optional
        Give the extension of the file you're looking for, default='.asc'    
    
    Returns
    -------
    selected_asciis : list
        List with strings of ASCII names
    """
    
    selected_asciis = []
    if str(glob.glob(location + '*' + 'Tr_01' + extension)).strip('[]') != "":  
        # if first case has Tr_01
        for i in range(nr_epochs):
            selected_asciis.append(str(glob.glob(location + '*' + 'Tr_' 
                    + "%02d" % (i+1) + extension)).strip('[]'))
    elif str(glob.glob(location + '*' + 'Tr_1' + extension)).strip('[]') != "":  
        # if first case has Tr_1
        for i in range(nr_epochs):
            selected_asciis.append(str(glob.glob(location + '*' 
                    + 'Tr_' + str(i+1) + extension)).strip('[]'))
    elif str(glob.glob(location + '*' + 'Tr01' + extension)).strip('[]') != "":  
        # if first case has Tr01
        for i in range(nr_epochs):
            selected_asciis.append(str(glob.glob(location + '*' + 'Tr' 
                    + "%02d" % (i+1) + extension)).strip('[]'))
    elif str(glob.glob(location + '*' + 'Tr1' + extension)).strip('[]') != "": 
        # if first case has Tr1
        for i in range(nr_epochs):
            selected_asciis.append(str(glob.glob(location + '*' + 'Tr' 
                    + str(i+1) + extension)).strip('[]'))
    return selected_asciis            


def cal_power_spectrum (timeseries, nr_rois=np.arange(92), fs=1250, 
            window='hamming', nperseg=4096, scaling='spectrum', 
            plot_figure=False, title_plot='average power spectrum'):

    """ Calculate (and plot) power spectrum of timeseries
    
    Parameters
    ----------
    timeseries: DataFrame with ndarrays
        Rows are timepoints, columns are rois/electrodes
    rois: int, optional
        Give list with rois/electrodes you want to include, 
        default=np.arange(92)
    fs: int, optional    
        Sample frequency, default=1250    
    window: str or tuple, optional
        Type of window you want to use, check spectral.py for details, 
        default='hamming'
    nperseg : int, optional    
        Length of each segment, default=4096
    scaling : str, optional
        'density' calculates the power spectral density (V**2/Hz), 'spectrum' 
        calculates the power spectrum (V**2), default='spectrum'
    plot_figure: bool
        Creates a figure of the mean + std over all rois/electrodes, 
        default=False
    title_plot: str
        Give title of the plot, default='average power spectrum'
        
    Returns
    -------
    f: ndarray
        Array with sample frequencies (x-asix of power spectrum plot)
    pxx: ndarray
        Columns of power spectra for each roi/VE
            
    """
  
    pxx = np.empty([int(nperseg/2+1), np.size(nr_rois)])
    
    i=0
    for roi in nr_rois:
        (f, pxx[:,i]) = signal.welch(timeseries[roi].values, fs, window, 
                                     nperseg, scaling=scaling)
        i=i+1
    if plot_figure==True:
        plt.figure()
        plt.plot(f, np.mean(pxx,1), color='teal')
        plt.plot(f, np.mean(pxx,1)+np.std(pxx,1), color='teal', linewidth=0.7)
        plt.plot(f, np.mean(pxx,1)-np.std(pxx,1), color='teal', linewidth=0.7)
        plt.fill_between(f, np.mean(pxx,1)+np.std(pxx,1), np.mean(pxx,1)
                         -np.std(pxx,1), color='teal', alpha=0.2)
        plt.xlim(0, 50)
        plt.xlabel('Frequency (Hz)')
        plt.title(title_plot)
        plt.show()
    return f, pxx


def find_nearest(array, value):
    """ Find nearest value of interest in array (used for frequencies, 
    no double value issues)
    
    Parameters
    ----------
    array: array
        Give the array in which you want to find index of value nearest-by
    value: int or float
        The value of interest
     
    Return
    ------
    idx: int
        Index of value nearest by value of interest    
    
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def run_loop_powerspectrum(subject_list, extension='.asc', nr_epochs=10, 
            nr_rois=np.arange(92), Fs=1250, window_length=4096, 
            freq_range=[0.5, 48]):

    """ Calculate power spectrum for all cases within the subject_list
    
    Parameters
    ----------
    subject_list: DataFrame, list, or string
        DataFrame should contain a column 'location_MEG' with the directory of 
        all MEG ASCIIs for each subject. List should contain only 1 column 
        describing the location. String should be the name of the ASCII
    extension: str, optional
        Give the extension of ASCIIs, '.txt' or '.asc', default='.asc', 
        only relevant when giving Dataframe or list
    nr_epochs: int, optional
        Give number of epochs you want to include, default=10,
        not relevant when input is single ascii
    nr_rois: int, optional
        Give list with rois you want to analyse, default=np.arange(92)
    Fs: int, optional
        Sample frequency, default=1250
    window_length: int, optional
        Window length to calculate power spectrum, default=4096
    freq_range: list, optional
        Gives the upper and lower boundaries to calculate the broadband power, 
        default=[0.5, 48]    
    
    Return
    ------
    mean_pxx: ndarray (size: len(subjects), len(power spectrum), nr_rois)
        Power spectum for all subjects, and all rois/VE, averaged over 
        nr_epochs
    broadband_power : ndarray (size: len(subjects), nr_rois)
        Broadband power between freq_range
    f: ndarray
        Array with sample frequencies (x-axis of power spectrum plot)
   
    """  
    
    mean_pxx = np.empty([len(subject_list), int(window_length/2+1), np.size(nr_rois)])
    broadband_power = np.empty([len(subject_list), np.size(nr_rois)])
    if isinstance(subject_list, str):
        single_ascii = subject_list
        timeseries =pd.read_csv(single_ascii, index_col=False, 
                                    header=None, delimiter='\t') 
        f, pxx = cal_power_spectrum(timeseries, nr_rois=nr_rois, fs=Fs, plot_figure=False, title_plot='power spectrum')
        mean_pxx= pxx            
        broadband_power = np.sum(mean_pxx[find_nearest(f,freq_range[0]):find_nearest(f, freq_range[1]),:], axis=0)
    else:
        if ~isinstance(subject_list, pd.DataFrame): 
            # check if subject_list is dataframe (or list)
            subject_list=pd.DataFrame(subject_list, columns = ['location_MEG']) 
            # save dataframe with column 'location'
            for sub in range(len(subject_list)):
                asciis = get_names_epochs(subject_list['location_MEG'][sub], nr_epochs, extension)
                sub_pxx = np.zeros((nr_epochs, int(window_length/2+1), np.size(nr_rois)))
                mean_pxx[sub,:,:] = 'nan'
                broadband_power[sub,:] = 'nan'
                f = np.empty(int(window_length/2+1))
                print(sub)
                if len(asciis) == 0:
                    print("no ASCIIs available for: " + subject_list['location_MEG'][sub])
                    continue    
                for j in range(nr_epochs):
                    timeseries =pd.read_csv(asciis[j][1:-1], index_col=False, 
                                    header=None, delimiter='\t') 
                    # First and last character should be removed otherwise 
                    # it's not working
                    # Compute power spectrum
                    f, pxx = cal_power_spectrum(timeseries, nr_rois=nr_rois, fs=Fs, 
                                                plot_figure=False, title_plot='average power spectrum - epoch: '+ str(j))
                    sub_pxx[j,:,:] = pxx           
                    mean_pxx[sub,:,:] = np.nanmean(sub_pxx, axis=0)    
                    broadband_power[sub,:] = np.sum(mean_pxx[sub,find_nearest(f,
                                    freq_range[0]):find_nearest(f,
                                    freq_range[1]),:], axis=0)
    return mean_pxx, broadband_power, f  


def cal_FOOOF_parameters(pxx, f, freq_range=[0.5, 48]):
    """ Obtain slope and offset using the FOOOF algorithm
    Reference paper: Haller M, Donoghue T, Peterson E, Varma P, Sebastian P, 
    Gao R, Noto T, Knight RT, Shestyuk A, Voytek B (2018) Parameterizing Neural
    Power Spectra. bioRxiv, 299859. doi: https://doi.org/10.1101/299859
    
    Reference Github: https://fooof-tools.github.io/fooof/index.html
    
    Parameters
    ----------
    pxx: ndarray
        Column of power spectra 
    f: 1d-array
        Array with sample frequencies (x-asix of power spectrum plot)
    freq_range: list, optional
        Gives the upper and lower boundaries to calculate the broadband power, 
        default=[0.5, 48]
    
    Returns
    -------
    FOOOF_offset: float
        Offset 
    FOOOF_slope: float
        Slope  
    
    """
    
    # initialize FOOOF oject
    fm = FOOOF()
    # create variables
    fm.fit(f, pxx, freq_range)
    FOOOF_offset = fm.background_params_[0]
    FOOOF_slope = fm.background_params_[1]
    time.sleep(1) # heavy algorithm
    
    return FOOOF_offset, FOOOF_slope


# -----------------------------------------------------------------------------

# set nice level to 10, especially FOOOF algorithm is heavy!
os.nice(10)


###########################
# Settings                #   
###########################

# 1. Select if you want to preprocess 1 single ASCII or multiple
process_single_ASCII='no' # either 'yes' or 'no'
# 2a. If you selected yes provide the location and name of the ASCII
# (e.g. '/data/KNW/NO-cohorten/Scans/sub-9110/meg/T1/AAL/9110_Tr01.asc' )
single_ascii = '/data/doorgeefluik/mumo_002_OD1_Tr01.asc'

# 2b_1. If you selected no; you can create an list of subjects you want to process
# an example is given here: '/data/KNW/t.numan/GOALS/example_MEG_list.csv'
# only available for NO-cohorten data!
subject_list = '/data/KNW/t.numan/GOALS/example_MEG_list.csv'
dir_input = '/data/KNW/NO-cohorten/Scans/'
# 2b_2. You can select AAL or BNA atlas
atlas = 'AAL' # AAL or BNA
# 2b_3. Select the number of epochs you want to preprocess
nr_epochs= 10 # number of epochs

# 3. Select which roi or rois you want to analyze
# if you want to analyze 1 roi, specify its number (nr_rois = (10,))
nr_rois = (10,) # only run for roi 11 # note that python indexes at 0!
# if you want to analyze multiple rois, create list with these rois
# (for example nr_rois = np.arange(78) for all 78 cortical AAL rois)

# 4. Set sample frequency (1250 Hz for Elekta data)
Fs = 1250 # sample frequency
# 5. Set frequency range you want to study
freq_range=[0.5, 48] # frequency range you want to analyze

# 6. Give output directory
dir_output = '/data/KNW/t.numan/GOALS/tmp/'
# 7a. Do you want to save the output? 
save_output = True # you can save output
# 7b. Provide relevant name to save your run
name_to_save = 'test_mumo_002_roi_11'


###########################
# Run analysis            # 
###########################

if process_single_ASCII=='no':
    # read the csv file containing Case_ID and MM of each included MEG
    # see example '/data/KNW/t.numan/GOALS/example_MEG_list.csv'
    subjects = pd.read_csv(subject_list, delimiter=';', header=0, usecols=[0,1])
    subjects['Case_ID'] = subjects['Case_ID'].apply(lambda x: str(x).zfill(4)) 
    # add padding zeros when necessary

    # get location of ASCIIs of MEG 
    subjects['location_MEG'] = dir_input + 'sub-'+ subjects['Case_ID'] + '/meg/' + subjects['MM'] + '/' + atlas + '/' 
    ###### Compute power spectrum ######
    ## takes a while when multiple subjects are analyzed at once
    mean_pxx, broadband_power, f = run_loop_powerspectrum(subjects, 
                    nr_epochs=nr_epochs, nr_rois=nr_rois, Fs=Fs, 
                    window_length=4096, freq_range=freq_range)
    
elif process_single_ASCII=='yes':
    mean_pxx, broadband_power, f = run_loop_powerspectrum(single_ascii, 
                    nr_epochs=1, nr_rois=nr_rois, Fs=Fs, 
                    window_length=4096, freq_range=freq_range)
else:
    print('wrong input process_single_ASCII')

    
# mean_pxx contains the average power spectra over nr_epochs for all subjects 
# for all rois, broadband_power is the area under the average power spectra 
# over the frequency range for all subjects and rois, f gives the frequencies 
# of the power spectrum (can be useful when plotting power spectrum)

# save output
if save_output == True:
    np.save(dir_output + name_to_save + '_mean_pxx', mean_pxx)
    np.save(dir_output + name_to_save + '_broadband_power', broadband_power)
    print('saving power spectra and broadband power')
    

###### Compute FOOOOF offset and slope ######
    
if process_single_ASCII=='no':    
    # create empty arrays to store offset and slope values
    offset = np.empty([mean_pxx.shape[0], mean_pxx.shape[2]]) 
    # index 0 -> number of subjects, index 2 -> number of rois
    slope = np.empty([mean_pxx.shape[0], mean_pxx.shape[2]])

    print('Running FOOOF.........................')

    # run across all subjects in your list
    for sub in range(len(subjects)):
        print('subject : ' + str(sub)) # print which subject is analyzed
        if np.isnan(mean_pxx[sub,1,0]): 
            # if there is no mean_pxx  for a subject, set offset/slope to nan
            offset[sub,] = np.nan
            slope[sub,] = np.nan
        else:
            i=0
            for roi in nr_rois:
                print('roi : ' + str(roi))
                offset[sub,i], slope[sub,i] = cal_FOOOF_parameters(mean_pxx[sub,:,i], 
                                                                       f, freq_range=[0.5, 48])
                i = i + 1
                time.sleep(0.05) 
                # by having a pause it should limit the %CPU continuously
                 # might be relevant when computing FOOOF for large dataset        
        if sub % 10 ==0 & save_output==True:
            print('save offset and slope in between after (another) 10 subjects')
            np.save(dir_output + name_to_save + '_FOOOF_offset', offset)
            np.save(dir_output + name_to_save + '_FOOOF_slope', slope)

 
elif process_single_ASCII=='yes':
    i=0
    offset = np.empty(mean_pxx.shape[1])
    slope = np.empty(mean_pxx.shape[1])
    for roi in nr_rois:
        print('roi : ' + str(roi))
        offset[i], slope[i] = cal_FOOOF_parameters(mean_pxx[:,i], f, freq_range=[0.5, 48])
        i = i + 1
        time.sleep(0.05) 
        # by having a pause it should limit the %CPU continuously
        # might be relevant when computing FOOOF for large dataset        
    
   
if save_output == True:          
    np.save(dir_output + name_to_save + '_FOOOF_offset', offset)
    np.save(dir_output + name_to_save + '_FOOOF_slope', slope)
    # save offset and slope output of FOOOF
    print("FOOOF-offset is saved: "+ dir_output + name_to_save 
          + "_FOOOF_offset")
    print("FOOOF-slope is saved: "+ dir_output + name_to_save + "_FOOOF_slope")

    
######## Combine results with Case_ID in 'subjects' variable ######
    # only whem multiple cases are processed at once 
if process_single_ASCII=='no':   
    subjects['mean_broadbandpower'] = np.mean(broadband_power, axis=1)     
    subjects['mean_offset'] = np.mean(offset, axis=1)
    subjects['mean_slope'] = np.mean(slope, axis=1)

    subjects.to_csv(dir_output + name_to_save + '_summary.csv' , sep=';', 
                index=False)
    # you can import this .csv file in SPSS or Excel
    

