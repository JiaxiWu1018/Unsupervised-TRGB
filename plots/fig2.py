#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 17:28:42 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from astropy.io import fits
from astropy.table import Table

def read_fits(filename):
    hdulist = fits.open(filename)
    dat = Table.read(hdulist[1], format = 'fits')
    photDF = dat.to_pandas()
    return photDF

plt.style.use('custom.mplstyle')
mpl.rcParams['text.usetex'] = False

#change version here
version = 4

info = pd.read_csv('../info/ghosts_info_v{:s}.csv'.format(str(version)))
judge = ['3031' in info['field'][i] for i in range(len(info))]
info = info[judge].reset_index(drop=True)
fig = plt.figure(figsize=(24, 30))

name = ['01','02','03','04','05','06','07','08','09','10','11']
for i in range(len(info)):
    ax = fig.add_subplot(3, 4, i+1)
    data = pd.read_csv('../csv/{:s}.csv'.format(info['field'][i]))
    color, mag = data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS']
    ax.plot(color, mag, 'k.', ms=2)
    ax.plot([-1, 3.5], [info['faint_lim'][i]]*2, 'g-', lw=4)
    ax.plot([info['red_lim'][i]]*2, [30, 20], 'r-', lw=4)
    
    m, b = info['slope_bc'][i], info['inter_bc'][i]
    yy = np.array([30, 20])
    xx1 = (yy - b) / m
    xx2 = (yy - b) / m + 1.0
    ax.plot(xx1, yy, c='orange', ls='--', lw=4)
    ax.plot(xx2, yy, c='orange', ls='--', lw=4)
    
    ax.text(1.65, 27.45, 'RA={:.2f}'.format(info['ra'][i]), fontsize=18)
    ax.text(1.65, 27.75, 'DEC={:.2f}'.format(info['dec'][i]), fontsize=18)
    ax.set_xlim(-1, 3.5)
    ax.set_ylim(28.25, 22.25)
    ax.set_title('Halo Field{:s}'.format(name[i]))
    
    if i in [0, 4, 8]:
        ax.set_ylabel('F814W')
    if i in [7, 8, 9, 10]:
        ax.set_xlabel('F606W-F814W')

#plt.savefig('CMD.png')
plt.show()