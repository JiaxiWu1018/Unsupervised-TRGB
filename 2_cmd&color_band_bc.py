#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 14:58:34 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def color_band(color, mag):
    num = 0
    for b in np.arange(25, 35, 0.1):
        for m in np.arange(-7, 0, 0.1):
            xx1 = (mag - b) / m
            xx2 = (mag - b) / m + 1.0
            judge = (color < xx2) & (color > xx1)
            if np.sum(judge) > num:
                num = np.sum(judge)
                slope, inter = m, b
    return slope, inter

info = pd.read_csv('ghosts_info_1.csv')
slope, inter = np.zeros(len(info)), np.zeros(len(info))
for i in range(len(info)):
    field = info['field'][i]
    print(field)
    data = pd.read_csv('csv/{:s}.csv'.format(field))
    color, mag = data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS']
    slope[i], inter[i] = color_band(color, mag)
info['slope_bc'], info['inter_bc'] = slope, inter

i = 0
while i < len(info):
    gal = info['galaxy'][i]
    judge = [(gal == info['galaxy'][j]) for j in range(len(info))]
    subinfo = info[judge].reset_index(drop=True)
    i += len(subinfo)
    
    fig, ax = plt.subplots(1, len(subinfo), figsize = (3.5*len(subinfo), 6))
    for j in range(len(subinfo)):
        data = pd.read_csv('csv/{:s}.csv'.format(subinfo['field'][j]))
        color, mag = data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS']
        red_lim = 0.3 + subinfo['A606'][j] - subinfo['A814'][j]
        faint_lim = subinfo['inter_bc'][j] + subinfo['slope_bc'][j] * red_lim
        ax[j].plot(color, mag, 'k.', ms=1)
        ax[j].plot([-1, 3.5], [faint_lim]*2, 'g-', lw=2)
        ax[j].plot([red_lim]*2, [35, 20], 'r-', lw=2)
        
        m, b = subinfo['slope_bc'][j], subinfo['inter_bc'][j]
        yy = np.array([35, 20])
        xx1, xx2 = (yy - b) / m, (yy - b) / m + 1.0
        ax[j].plot(xx1, yy, c='orange', ls='--', lw=2)
        ax[j].plot(xx2, yy, c='orange', ls='--', lw=2)
        
        ax[j].set_xlim(-1, 3.5)
        blim = np.round(np.max(mag)*2)/2
        ax[j].set_ylim(blim, blim-6)
        ax[j].set_title(subinfo['field'][j])
    plt.savefig('cmd/{:s}.png'.format(gal))
    plt.show()

info.to_csv('ghosts_info_2.csv', index=False)
