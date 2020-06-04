#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 15:04:02 2020

@author: t.numan
"""
# reviewed by Linda Douw 20200603
# updated by Tianne Numan 20200604

import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import signal
import pandas as pd
from fooof import FOOOF

def cal_power_spectrum (timeseries, rois, fs, window='hamming', nperseg=4096, scaling='spectrum', plot_figure=False, title_plot='average power spectrum'):
    pxx = np.empty([int(nperseg/2+1), len(rois)])
    for roi in rois:
        (f, pxx[:,roi]) = scipy.signal.welch(timeseries[roi].values, fs, window, nperseg, scaling=scaling)
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

# you can enter the output of cal_power_spectrum, pxx is matrix with a power spectrum for all rois, f is the vector with frequencies
# rois = number of rois, default=210
# freq_range = the frequency range you want to fit the model, default=[0.5, 48]

# also check FOOOF algorithm: fooof-tools.github.io/fooof/auto_tutorials/index.html

def cal_FOOOF_parameters(pxx, f, rois=210, freq_range=[0.5, 48]):
    # initialize FOOOF oject
    fm = FOOOF()
    # create variables
    FOOOF_offset = np.empty(rois)
    FOOOF_slope = np.empty(rois)
    # loop over all rois 
    for roi in list(range(rois)):
            # model the power spectrum with FOOOF for each roi in rois
            # LD: this part does not run because the size of the freqs var is somehow not consistent
            fm.fit(f, pxx[:,roi], freq_range)
            FOOOF_offset[roi] = fm.background_params_[0]
            FOOOF_slope[roi] = fm.background_params_[1]
    
    return FOOOF_offset, FOOOF_slope


# Example
dir_MEG='/data/KNW/NO-cohorten/Scans/sub-0015/meg/T1/BNA/'   
name_ASCII='1_100_WITH_200_WITH_246_VE_12.000to25.106_Tr_1.asc'

timeseries = pd.read_csv(dir_MEG+name_ASCII, index_col=False, header=None, delimiter='\t')
f, pxx = cal_power_spectrum(timeseries,list(range(210)), 1250, plot_figure=True)

# enter power spectrum in FOOOF algorithm, run for the first 1 only
offset, slope = cal_FOOOF_parameters(pxx, f,rois=1)

