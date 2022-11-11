#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 18:30:42 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

def count_field(det):
    field, gal = [], []
    for i in range(len(det)):
        if det['field'][i] not in field:
            field.append(det['field'][i])
        if det['galaxy'][i] not in gal:
            gal.append(det['galaxy'][i])
    return field, gal

plt.style.use('custom.mplstyle')
mpl.rcParams['text.usetex'] = False

#change version here
version = 4.1

color_list = ['#7e1e9c', '#15b01a', '#0343df', '#ff81c0', '#653700', '#e50000', '#95d0fc', '#029386', '#f97306', '#96f97b', '#c20078', '#ffff14',
              '#75bbfd', '#929591', '#89fe05', '#bf77f6', '#9a0eea', '#033500', '#06c2ac', '#c79fef', '#00035b', '#d1b26f', '#00ffff', '#13eac9',
              '#06470c', '#ae7181', '#35063e', '#01ff07', '#650021', '#6e750e', '#ff796c', '#e6daa6', '#0504aa', '#001146', '#cea2fd', '#000000',
              '#ff028d', '#ad8150', '#c7fdb5', '#ffb07c', '#677a04', '#cb416b', '#8e82fe', '#53fca1', '#aaff32', '#380282', '#ceb301', '#ffd1df',
              '#cf6275', '#0165fc', '#0cff0c', '#c04e01', '#04d8b2', '#01153e', '#3f9b0b', '#d0fefe', '#840000', '#be03fd', '#c0fb2d', '#a2cffe',
              '#dbb40c', '#8fff9f', '#580f41', '#4b006e', '#8f1402', '#014d4e', '#610023', '#aaa662', '#137e6d', '#7af9ab', '#02ab2e', '#9aae07']
fn_dict = {'ngc0253':4,'ngc0891':4,'ngc2403':5,'ngc3031':11,'ngc4244':4,
            'ngc4565':3,'ngc4736':2,'ngc5236':10,'ngc7793':4,'ngc7814':3}

det = pd.read_csv('../detection/ghosts_detection_v{:s}.csv'.format(str(version)))

fig, ax = plt.subplots(figsize = (20, 12))
ax.plot([-1, 10], [0]*2, 'k-', lw=4)
i, count, gal = 0, 0, []
while i < len(det):
    judge = [det['galaxy'][j] == det['galaxy'][i] for j in range(len(det))]
    subdet = det[judge].reset_index(drop=True)
    if len(subdet) <= 1:
        i += len(subdet)
        continue
    i += len(subdet)
    
    judge = (subdet['RGB AGB Ratio']>=4) & (subdet['# star below tip']>=200)
    subdet = subdet[judge].reset_index(drop=True)
    if len(subdet) <= 1:
        continue
    mean = np.mean(np.array(subdet['TRGB']))
    res = np.array(subdet['TRGB']) - mean
    gal.append(det['galaxy'][i-len(subdet)])
    ax.plot([count+0.5]*len(subdet), res, c=color_list[count], marker = 'D', ms=16, ls='')
    
    field, galname = count_field(subdet)
    per = len(field) / fn_dict[galname[0]] * 100
    tpf = len(subdet) / len(field)
    ax.text(count+0.2, -0.14, '${:.1f}\%$'.format(per), fontsize=24, c='red')
    ax.text(count+0.3, -0.18, '${:.1f}$'.format(tpf), fontsize=24, c='red')
    count += 1

ax.text(-0.45, -0.14, '$P_{valid}$', fontsize=24, c='red')
ax.text(-0.45, -0.18, '$N_{TPF}$', fontsize=24, c='red')

for i in range(len(gal)):
    gal[i] = 'NGC' + gal[i][3:]

ax.set_xlim(-0.5, len(gal))
ax.set_xticks(np.arange(0.5, len(gal), 1))
ax.set_xticklabels(gal, rotation=30)
ax.set_ylabel('Residual to Mean (mag)')
ax.set_ylim(-0.22, 0.25)
ax.set_yticks(np.arange(-0.20, 0.25, 0.1))
ax.grid()

#plt.savefig('variance.png')
plt.show()
