#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 15:04:02 2020

@author: t.numan
"""
# reviewed by Linda Douw 20200603

import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import signal
import pandas as pd
from fooof import FOOOF

# you can enter the output of cal_power_spectrum, pxx is matrix with a power spectrum for all rois, f is the vector with frequencies
# rois = list of rois (e.g. range(92) or a list with rois-names)
# freq_range = the frequency range you want to fit the model

# also check FOOOF algorithm: fooof-tools.github.io/fooof/auto_tutorials/index.html

# initialize FOOOF oject
fm = FOOOF()

def cal_FOOOF_parameters(pxx, f, rois, freq_range=[0.5, 48]):
    # create variables
    FOOOF_offset = np.empty([len(rois)])
    FOOOF_slope = np.empty([len(rois)])
    # loop over all rois 
    for roi in rois:
            # model the power spectrum with FOOOF for each roi in rois
            # LD: this part does not run because the size of the freqs var is somehow not consistent
            fm.fit(f, pxx[roi,], freq_range)
            FOOOF_offset[roi] = fm.background_params_[0]
            FOOOF_slope[roi] = fm.background_params_[1]
    
    return FOOOF_offset, FOOOF_slope