# function to calculate and plot power spectrum
# give at least timeseries, number of rois (rois) and sample frequency (fs) as input in the function.
# other parameters are optional
# Code review by Linda Douw on 29052020

# packages needed (?not sure if they are all needed for this script)
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import signal
import pandas as pd

# we read in the timeries here first, perhaps we need to do this in a wrapper script when using this function per project? 
timeseries=pd.read_csv('/data/doorgeefluik/mumo_002_OD1_Tr01.asc', index_col=False, header=None, delimiter='\t')

def cal_power_spectrum (timeseries, rois, fs, window='hamming', nperseg=4096, scaling='spectrum', plot_figure=False, title_plot='average power spectrum'):
    pxx = np.empty([int(nperseg/2+1), rois])
    #for roi in rois:
    for roi in range(rois):
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
