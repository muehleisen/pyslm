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
Ka = (2*np.pi*fa4)**2 * (10**(A1000/20))
Za = [0, 0, 0, 0]
Pa = 2*np.pi*np.array([fa4, fa4, fa3, fa2, fa1, fa1])


# C weighting filter has poles at 20.6 and 12194 Hz
C1000 = 0.0619
Kc = (2*np.pi*fa4)**2 * (10**(C1000/20))
Zc = [0, 0]
Pc = 2*np.pi*np.array([fa4, fa1, fa4, fa1])

FS_list = [22050, 44100, 48000, 96000]

# tolerance for Type 0 and Type 1 filters from S1.42
ftol = [10, 12.5, 16, 20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250,
        315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000,
        5000, 6300, 8000, 10000, 12500, 16000, 20000]
up1 = np.array([3, 2.5, 2, 2, 2, 1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 0.7, 1, 1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 2, 2, 2.5, 3])
low1 = np.array([np.inf, np.inf, 4, 2, 1.5, 1.5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                1, 1, 1, 1, 0.7, 1, 1, 1, 1, 1, 1, 1.5, 2, 2.5, 3, 5, 8, np.inf])
up0 = np.array([2, 2, 2, 2, 1.5, 1, 1, 1, 1, 1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7,
                0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 1, 1, 1, 2, 2, 2, 2])
low0 = np.array([5, 3, 3, 2, 1.5, 1, 1, 1, 1, 1, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7,
                0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 1, 1.5, 2, 3, 3, 3, 3])

# FS_list = [4]
for fs in FS_list:
    f = np.logspace(1, np.log10(30000), 1024)
    wn = f * 2 * np.pi/fs
    [ba, aa] = vf.awt_design(fs)
    [bc, ac] = vf.cwt_design(fs)
    # compute the frequency response of the analog weighting filter
    [wa, ha] = sg.freqs_zpk(Za, Pa, Ka, wn*fs)
    [wc, hc] = sg.freqs_zpk(Zc, Pc, Kc, wn*fs)

    [wza, hza] = sg.freqz(ba, aa, wn)
    [wzc, hzc] = sg.freqz(bc, ac, wn)
# end for fs in FS_list

    hadb = 20*np.log10(np.abs(ha))
    hcdb = 20*np.log10(np.abs(hc))
    hzadb = 20*np.log10(np.abs(hza))
    hzcdb = 20*np.log10(np.abs(hzc))

    diffa = hzadb - hadb
    diffc = hzcdb - hcdb

    plt.figure(figsize=[10, 6])

    plt.subplot(2, 2, 1)
    plt.semilogx(f, hzadb, f, hadb)

    plt.axis([10, 30000, -70, 10])
    plt.title(' A-weighting Filter for Fs = {}'.format(fs))

    plt.subplot(2, 2, 3)

    plt.semilogx(f, diffa, ftol, up0, 'r--', ftol, -low0, 'r--')
    plt.axis([10, 30000, -5, 5])

    plt.subplot(2, 2, 2)
    plt.semilogx(f, hzcdb, f, hcdb)

    plt.axis([10, 30000, -70, 10])
    plt.title(' C-weighting Filter for Fs = {}'.format(fs))

    plt.subplot(2, 2, 4)

    plt.semilogx(f, diffc, ftol, up0, 'r--', ftol, -low0, 'r--')
    plt.axis([10, 30000, -5, 5])

    plt.show()
