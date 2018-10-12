# test_SLM_filters.py
# This file tests the digital IIR filters defined VSLM_filter.py
#
# It compares the IIR filter to the analog filter response and plots the
# two responses with the type 1 filter tolerances
#
#

import scipy.signal as sg
import numpy as np
import matplotlib.pyplot as plt

import VSLM_filters as vf

# define poles of A weighted filter as per ANSI S1.42 annex
fa1 = 20.598997
fa2 = 107.65265
fa3 = 737.86223
fa4 = 12194.217

A1000 = 1.9997
K = (2*np.pi*fa4)**2 * (10**(A1000/20))
Z = [0, 0, 0, 0]
P = 2*np.pi*np.array([fa4, fa4, fa3, fa2, fa1, fa1])
[Ba, Aa] = sg.zpk2tf(Z, P, K)

FS_list = [22050, 44100, 48000, 96000]




# FS_list = [4]
for fs in FS_list:
    f = np.logspace(1, np.log10(fs/2), 512)
    wn = f * 2 * np.pi/fs
    [b, a] = vf.awt_design(fs)

    # compute the frequency response of the actual a weighted filter
    [wa, ha] = sg.freqs_zpk(Z, P, K, wn*fs )
    [wz, hz] = sg.freqz(b, a, wn)
# end for fs in FS_list

# print(ha)
    hadb = 20*np.log10(np.abs(ha))
    hzdb = 20*np.log10(np.abs(hz))

    plt.semilogx(f, hzdb, f, hadb)
    plt.show()
# print(ha)
# print(hz)

