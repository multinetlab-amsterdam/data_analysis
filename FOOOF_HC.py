# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:25:04 2019

@author: Tianne
"""

# use FOOOF to calculate model with aperiodic part to calculate power spectrum
# import power-spectrum calculated with Matlab

import numpy as np
# import mne
import matplotlib.pyplot as plt
#import tqdm
import scipy
from scipy import stats
import sys
from scipy.integrate import simps
from fooof import FOOOF
import os, glob


# define frequency range across which to model the spectrum
freq_range = [0.5, 48]
bna_rois = list(range(210))
Fs = 1250
#freqs_origial = np.arange(0,Fs/2, 625/2049)

# initialize FOOOF oject
fm = FOOOF()

dir_files = '/data/KNW/t.numan/GOALS/retrospective_part/matlab_output/HC_spectra/'


# import frequency axis
freqs_file = 'MEG_f_powerspectrum.txt'
freqs = np.loadtxt(dir_files+freqs_file)


# get all .txt files and save in sub_list
all_files = os.listdir(dir_files)
sub_list = []
for names in all_files:
    if"mean_powerspectra_Case" in names:
        sub_list.append(names)

# define empty variables for original power and flat power spectrum
AUC_abs_flat_spectrum = np.zeros([len(sub_list), 210])
AUC_abs_orig_spectrum = np.zeros([len(sub_list), 210])
flat_spectrum = np.zeros([210, 156])
orig_spectrum = np.zeros([210, 156])
one_over_f_spectrum = np.zeros([210, 156])
AUC_matlab_spectrum = np.zeros([len(sub_list), 210])
mean_flat_spectrum = np.zeros([len(sub_list), 156])
std_flat_spectrum = np.zeros([len(sub_list), 156])
mean_std1_flat_spectrum = np.zeros([len(sub_list), 156])
background_params_offset = np.zeros([len(sub_list), 210])
background_params_slope = np.zeros([len(sub_list), 210])


mean_orig_spectrum = np.zeros([len(sub_list), 156])
std_orig_spectrum = np.zeros([len(sub_list), 156])
mean_one_over_f = np.zeros([len(sub_list), 156])
std_one_over_f = np.zeros([len(sub_list), 156])


sub = 0
for file in sub_list:
    spectrum = np.loadtxt(dir_files+file)
    for roi in bna_rois:
        # model the power spectrum with FOOOF for each roi
       fm.fit(freqs, spectrum[roi, :], freq_range)
       flat_spectrum[roi, ] = fm._spectrum_flat
       orig_spectrum[roi, ] = fm.power_spectrum
       one_over_f_spectrum[roi, ] = fm._bg_fit
       AUC_abs_flat_spectrum[sub, roi] = sum(abs(fm._spectrum_flat))
       AUC_abs_orig_spectrum[sub, roi] = sum(abs(fm.power_spectrum))
       AUC_matlab_spectrum[sub, roi] = sum(spectrum[roi, 2:157])
       # get offset and slope for each roi and subject
       background_params_offset[sub, roi] = fm.background_params_[0]
       background_params_slope[sub, roi] = fm.background_params_[1]
    
    mean_flat_spectrum[sub, ] = np.mean(flat_spectrum, axis=0)
    std_flat_spectrum[sub, ] = np.std(flat_spectrum, axis=0)
    mean_std1_flat_spectrum[sub, ] = mean_flat_spectrum[sub, ] + std_flat_spectrum[sub, ]

#    plt.figure()
#    plt.plot(fm.freqs, mean_flat_spectrum[sub, ], label=sub, color='k')
#    plt.plot(fm.freqs, mean_std1_flat_spectrum[sub, ], label='std', color='k')
#    plt.plot(fm.freqs, mean_flat_spectrum[sub, ] - std_flat_spectrum[sub, ], label='std', color='k')
#    plt.fill_between(fm.freqs, mean_std1_flat_spectrum[sub, ], mean_flat_spectrum[sub, ] - std_flat_spectrum[sub, ], facecolor='dodgerblue', alpha = 0.5)
#    plt.xlabel('frequency (Hz)')
#    plt.title('Average FOOOFed spectrum of subject ' + str(sub))
    
    mean_one_over_f[sub, ] = np.mean(one_over_f_spectrum, axis=0)
    std_one_over_f[sub, ] = np.std(one_over_f_spectrum, axis=0)
    
    
#    plt.figure()
#    plt.plot(fm.freqs, mean_one_over_f[sub, ], label=sub, color='r')
#    plt.plot(fm.freqs, mean_one_over_f[sub, ] + std_one_over_f[sub, ], label='std', color='r')
#    plt.plot(fm.freqs, mean_one_over_f[sub, ] - std_one_over_f[sub, ], label='std', color='r')
#    plt.fill_between(fm.freqs, mean_one_over_f[sub, ] + std_one_over_f[sub, ], mean_one_over_f[sub, ] - std_one_over_f[sub, ], facecolor='lightcoral', alpha=0.5)
#    plt.xlabel('frequency (Hz)')
#    plt.title('1/f plot of subject ' + str(sub))
    print(sub) 
    
    mean_orig_spectrum[sub, ] = np.mean(orig_spectrum, axis=0)
    std_orig_spectrum[sub, ] = np.std(orig_spectrum, axis=0)
    
#    plt.figure()
#    plt.plot(fm.freqs, mean_orig_spectrum[sub, ], label=sub, color='k')
#    plt.plot(fm.freqs, mean_orig_spectrum[sub, ] + std_orig_spectrum[sub, ], label='std', color='k')
#    plt.plot(fm.freqs, mean_orig_spectrum[sub, ] - std_orig_spectrum[sub, ], label='std', color='k')
#    plt.fill_between(fm.freqs, mean_orig_spectrum[sub, ] + std_orig_spectrum[sub, ], mean_orig_spectrum[sub, ] - std_orig_spectrum[sub, ], facecolor='forestgreen', alpha = 0.5)
#    plt.xlabel('frequency (Hz)')
#    plt.title('Original frequency spectrum of subject ' + str(sub))
    
    sub= sub + 1

# np.corrcoef(AUC_abs_flat_spectrum, AUC_abs_original)
corr_orig_FOOOFed = np.zeros([len(sub_list), 2])


from scipy.stats.stats import pearsonr

sub = 0
for file in sub_list:
    corr_orig_FOOOFed[sub, ] = pearsonr(AUC_abs_flat_spectrum[sub, ], AUC_abs_orig_spectrum[sub, ])
   # corr_orig_FOOOFed = tmp[0]
   # corr_p_orig_FOOOFed.append = tmp[1]
    plt.figure()
    plt.scatter(AUC_abs_flat_spectrum[sub, ], AUC_abs_orig_spectrum[sub, ], label=sub)
    plt.xlabel('FOOOFed power')
    plt.ylabel('Original power')
    plt.title('Correlation FOOOFed and original power of HC')
    plt.legend()
    
    plt.figure()
    plt.subplot(1,2,1)
    plt.scatter(AUC_matlab_spectrum[sub, ], AUC_abs_orig_spectrum[sub, ], label=sub)
    plt.xlabel('Matlab broadband power')
    plt.ylabel('Original broadband power')
    plt.subplot(1,2,2)
    plt.scatter(AUC_matlab_spectrum[sub, ], AUC_abs_flat_spectrum[sub, ], label=sub)
    plt.xlabel('Matlab broadband power')
    plt.ylabel('FOOOFed power')
    
    sub = sub + 1
    
     
 # create scatterplot of background parameters over the 11 subjects and 210 rois
 # repeat the bna rois 11 times, so that you can scatter the same amount of datapoints for x and y
# x=np.repeat(bna_rois,background_params_offset.shape[0])
# plt.figure()
# plt.scatter(x, background_params_offset, )

# BACKGROUND PARAMETERS - OFFSET
 mean_offset = np.mean(background_params_offset,0)
 mean_offset_sub= np.mean()
plt.figure()

sub =0
for file in sub_list:
    plt.figure(1)
    plt.scatter(bna_rois, background_params_offset[sub, ], label=sub)
    plt.xlabel('bna rois')
    plt.ylabel('Offset')
    plt.title('Background parameter - offset')
    plt.figure(2)
    plt.plot(bna_rois, background_params_offset[sub, ], label=sub)
       
    sub = sub + 1
    plt.figure(1)
plt.scatter(bna_rois, mean_offset, label='mean', s=40, facecolor='k')
#plt.legend(loc=7)

# BACKGROUND PARAMETERS - SLOPE
mean_slope = -np.mean(background_params_slope,0)
plt.figure()
sub =0
for file in sub_list:
    plt.scatter(bna_rois, -background_params_slope[sub, ], label=sub)
    plt.xlabel('bna rois')
    plt.ylabel('Slope')
    plt.title('Background parameter - slope')
       
    sub = sub + 1
plt.scatter(bna_rois,mean_slope, label='mean', s=40, facecolor='k')

