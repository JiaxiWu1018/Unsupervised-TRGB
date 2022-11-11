#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 12:29:13 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

def c4(n):
    k = n//2
    if n%2 == 0:
        return np.sqrt(2/(np.pi*(2*k-1)))*(2**(2*k-2))*(math.factorial(k-1)**2)/math.factorial(2*k-2)
    else:
        return np.sqrt(np.pi/k)*math.factorial(2*k-1)/(2**(2*k-1))/(math.factorial(k-1)**2)

plt.style.use('custom.mplstyle')
mpl.rcParams['text.usetex'] = False

#change version here
version, version_det = 4, 4.1

data = pd.read_csv('../detection/ghosts_detection_v{:s}.csv'.format(str(version_det)))
data['TRGB'] = data['TRGB'] - data['Extinction'] - data['A814']
info = pd.read_csv('../info/ghosts_info_v{:s}.csv'.format(str(version)))
gallist = ['0253','0891','2403','3031','4244','4565','4736','5236','7793','7814']

nbtcut = np.array([50, 100, 200, 400])
racut = np.arange(1.5, 7, 0.5)
fig, ax = plt.subplots(3, 2, figsize=(36, 18))
line = []
for i in range(len(nbtcut)):
    data = data[data['# star below tip']>=nbtcut[i]].reset_index(drop=True)
    stdra, per, tpf, meanra, meannbt = [], [], [], [], []
    for j in range(len(racut)):
        datara = data[(data['RGB AGB Ratio']>=racut[j])&(data['RGB AGB Ratio']<=racut[j]+3)].reset_index(drop=True)
        meanra.append(np.mean(datara['RGB AGB Ratio']))
        meannbt.append(np.mean(datara['# star below tip']))
        if (i == 0) & (j == 0):
            field = []
            for k in range(len(datara)):
                if datara['field'][k] not in field:
                    field.append(datara['field'][k])
            numoffield = len(field)
        std, num, field = [], [], []
        for k in range(len(datara)):
            if datara['field'][k] not in field:
                field.append(datara['field'][k])
        for k in range(len(gallist)):
            gal = gallist[k]
            judge = [gal in datara['field'][m] for m in range(len(datara))]
            if np.sum(judge) >= 2:
                std.append(np.std(datara['TRGB'][judge])/c4(int(np.sum(judge))))
                num.append(int(np.sum(judge)))
        summ = 0
        for k in range(len(std)):
            summ += std[k] * num[k]
        stdra.append(summ / np.sum(np.array(num)))
        per.append(100 * len(field) / numoffield)
        tpf.append(len(datara)/len(field))
    
    ax[0,0].plot(racut, stdra, marker='*', ms=15, lw=4, label = '$N \geq {:d}$'.format(nbtcut[i]))
    ax[1,0].plot(racut, per, marker='*', ms=15, lw=4)
    ax[2,0].plot(racut, tpf, marker='*', ms=15, lw=4)
    ret, = ax[0,1].plot(meanra, stdra, marker='*', ms=15, lw=4, label = '$N \geq {:d}$'.format(nbtcut[i]))
    line.append(ret)
    ax[1,1].plot(meanra, per, marker='*', ms=15, lw=4)
    ax[2,1].plot(meanra, tpf, marker='*', ms=15, lw=4)
    
err = []
for j in range(len(meanra)):
    temp = np.exp(1.5*(3-meanra[j])) / (np.exp(1.5*(3-meanra[j]))+1) * (1/(meannbt[j]-100))**0.1 * 2
    temp = temp ** 2 + 0.04 ** 2
    err.append(np.sqrt(temp))
fit, = ax[0,1].plot(meanra, err, c='purple', ls='-.', lw=2, label = 'fit of Eq. (11)')

ax[0,0].set_xlim(1.2, 7)
ax[1,0].set_xlim(1.2, 7)
ax[2,0].set_xlim(1.2, 7)

ax[0,0].xaxis.set_major_locator(MultipleLocator(1))
ax[0,0].xaxis.set_minor_locator(MultipleLocator(0.5))
ax[0,0].yaxis.set_major_locator(MultipleLocator(0.2))
ax[0,0].yaxis.set_minor_locator(MultipleLocator(0.05))
ax[0,0].set_ylim(0.01, 1)
ax[0,0].set_yscale('log')
ax[0,0].set_yticks([1e-2,5e-2,1e-1,5e-1,1])
ax[0,0].set_yticklabels(['0.01','0.05', '0.10','0.50','1.00'])
ax[0,0].plot([1.2, 7], [0.1]*2, c='grey', ls='--', lw=2)
ax[0,0].plot([1.2, 7], [0.05]*2, c='grey', ls='--', lw=2)
ax[0,0].legend(loc = 'upper right')
ax[0,0].set_title('Multifield TRGB Dispersion vs Minimum Ratio')
ax[0,0].set_ylabel('$\sigma_{tot}$ (mag)')

ax[1,0].xaxis.set_major_locator(MultipleLocator(1))
ax[1,0].xaxis.set_minor_locator(MultipleLocator(0.5))
ax[1,0].yaxis.set_major_locator(MultipleLocator(20))
ax[1,0].yaxis.set_minor_locator(MultipleLocator(5))
ax[1,0].set_ylim(20, 105)
ax[1,0].plot([1.2, 7], [50]*2, c='grey', ls='--', lw=2)
ax[1,0].plot([1.2, 7], [30]*2, c='grey', ls='--', lw=2)
ax[1,0].set_ylabel('$P_{valid}$')

ax[2,0].xaxis.set_major_locator(MultipleLocator(1))
ax[2,0].xaxis.set_minor_locator(MultipleLocator(0.5))
ax[2,0].yaxis.set_major_locator(MultipleLocator(0.5))
ax[2,0].yaxis.set_minor_locator(MultipleLocator(0.1))
ax[2,0].set_ylim(0.9, 2.5)
ax[2,0].plot([1.2, 7], [1.5]*2, c='grey', ls='--', lw=2)
ax[2,0].plot([1.2, 7], [1]*2, c='grey', ls='--', lw=2)
ax[2,0].set_ylabel('$N_{TPF}$')
ax[2,0].set_xlabel('Minimum $R$')

ax[0,1].set_xlim(2, 8)
ax[1,1].set_xlim(2, 8)
ax[2,1].set_xlim(2, 8)

ax[0,1].xaxis.set_major_locator(MultipleLocator(1))
ax[0,1].xaxis.set_minor_locator(MultipleLocator(0.5))
ax[0,1].yaxis.set_major_locator(MultipleLocator(0.2))
ax[0,1].yaxis.set_minor_locator(MultipleLocator(0.05))
ax[0,1].set_ylim(0.01, 1)
ax[0,1].set_yscale('log')
ax[0,1].set_yticks([1e-2,5e-2,1e-1,5e-1,1])
ax[0,1].set_yticklabels(['0.01','0.05', '0.10','0.50','1.00'])
ax[0,1].plot([2, 8], [0.1]*2, c='grey', ls='--', lw=2)
ax[0,1].plot([2, 8], [0.05]*2, c='grey', ls='--', lw=2)
legend1 = ax[0,1].legend(handles = line, loc = 'upper right')
ax[0,1].add_artist(legend1)
ax[0,1].legend(handles = [fit], loc = 'lower left')
ax[0,1].set_title('Multifield TRGB Dispersion vs Mean Ratio')

ax[1,1].xaxis.set_major_locator(MultipleLocator(1))
ax[1,1].xaxis.set_minor_locator(MultipleLocator(0.5))
ax[1,1].yaxis.set_major_locator(MultipleLocator(20))
ax[1,1].yaxis.set_minor_locator(MultipleLocator(5))
ax[1,1].set_ylim(20, 105)
ax[1,1].plot([2, 8], [50]*2, c='grey', ls='--', lw=2)
ax[1,1].plot([2, 8], [30]*2, c='grey', ls='--', lw=2)

ax[2,1].xaxis.set_major_locator(MultipleLocator(1))
ax[2,1].xaxis.set_minor_locator(MultipleLocator(0.5))
ax[2,1].yaxis.set_major_locator(MultipleLocator(0.5))
ax[2,1].yaxis.set_minor_locator(MultipleLocator(0.1))
ax[2,1].set_ylim(0.9, 2.5)
ax[2,1].plot([2, 8], [1.5]*2, c='grey', ls='--', lw=2)
ax[2,1].plot([2, 8], [1]*2, c='grey', ls='--', lw=2)
ax[2,1].set_xlabel('Mean $R$')

#plt.savefig('Dispersion.png')
plt.show()
