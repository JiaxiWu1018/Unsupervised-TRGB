#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 18:25:23 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib as mpl

def c4(n):
    k = n//2
    if n%2 == 0:
        return np.sqrt(2/(np.pi*(2*k-1)))*(2**(2*k-2))*(math.factorial(k-1)**2)/math.factorial(2*k-2)
    else:
        return np.sqrt(np.pi/k)*math.factorial(2*k-1)/(2**(2*k-1))/(math.factorial(k-1)**2)

plt.style.use('custom.mplstyle')
mpl.rcParams['text.usetex'] = False

fig, ax = plt.subplots(1, 2, figsize=(32, 12))

det = pd.read_csv('../detection/ghosts_detection_v1.2.csv')
judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
TRGB = np.array(det['TRGB'])
res = TRGB - np.mean(TRGB)
std = np.std(res)/c4(len(res))
ax[0].plot([1]*len(res), res, ls='', marker='D', ms=16, c='orange')
ax[0].plot([0.5, 1.5], [std]*2, 'g--', lw=3)
ax[0].plot([0.5, 1.5], [-std]*2, 'g--', lw=3)

det = pd.read_csv('../detection/ghosts_detection_v4.3.csv')
judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
TRGB = np.array(det['TRGB'])
res = TRGB - np.mean(TRGB)
std = np.std(res)/c4(len(res))
ax[0].plot([3]*len(res), res, ls='', marker='D', ms=16, c='orange')
ax[0].plot([2.5, 3.5], [std]*2, 'g--', lw=3)
ax[0].plot([2.5, 3.5], [-std]*2, 'g--', lw=3)

det = pd.read_csv('../detection/ghosts_detection_v4.4.csv')
judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
TRGB = np.array(det['TRGB'])
res = TRGB - np.mean(TRGB)
std = np.std(res)/c4(len(res))
ax[0].plot([5]*len(res), res, ls='', marker='D', ms=16, c='orange')
ax[0].plot([4.5, 5.5], [std]*2, 'g--', lw=3)
ax[0].plot([4.5, 5.5], [-std]*2, 'g--', lw=3)

judge = det['RGB AGB Ratio'] >= 4
det = det[judge].reset_index(drop=True)
TRGB = np.array(det['TRGB'])
res = TRGB - np.mean(TRGB)
std = np.std(res)/c4(len(res))
ax[0].plot([7]*len(res), res, ls='', marker='D', ms=16, c='orange')
ax[0].plot([6.5, 7.5], [std]*2, 'g--', lw=3)
ax[0].plot([6.5, 7.5], [-std]*2, 'g--', lw=3)

judge = det['# star below tip'] >= 200
det = det[judge].reset_index(drop=True)
TRGB = np.array(det['TRGB'])
res = TRGB - np.mean(TRGB)
std = np.std(res)/c4(len(res))
ax[0].plot([9]*len(res), res, ls='', marker='D', ms=16, c='orange')
ax[0].plot([8.5, 9.5], [std]*2, 'g--', lw=3)
ax[0].plot([8.5, 9.5], [-std]*2, 'g--', lw=3)

det = pd.read_csv('../detection/ghosts_detection_v4.1.csv')
judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
judge = (det['RGB AGB Ratio'] >= 4) & (det['# star below tip'] >= 200)
det = det[judge].reset_index(drop=True)
TRGB = np.array(det['TRGB'])
res = TRGB - np.mean(TRGB)
std = np.std(res)/c4(len(res))
ax[0].plot([11]*len(res), res, ls='', marker='D', ms=16, c='orange')
ax[0].plot([10.5, 11.5], [std]*2, 'g--', lw=3)
ax[0].plot([10.5, 11.5], [-std]*2, 'g--', lw=3)

ax[0].set_xticks(np.arange(1, 12, 2))
ax[0].set_xticklabels(['raw', 'spatial', 'color band', '$R$ cut', '$N_{-,1.0}$ cut', 'smoothing'],fontsize=22)
ax[0].set_ylabel('Residual to Mean (mag)')
ax[0].set_ylim(-3.5, 3.5)
ax[0].set_yticks(np.arange(-3, 4, 1))
ax[0].grid()

det = pd.read_csv('../detection/ghosts_detection_v1.2.csv')
judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
TRGB = np.array(det['TRGB'])
res = TRGB - np.mean(TRGB)
std = np.std(res)/c4(len(res))
ax[1].plot([1]*len(res), res, ls='', marker='D', ms=16, c='orange')
ax[1].plot([0.5, 1.5], [std]*2, 'g--', lw=3)
ax[1].plot([0.5, 1.5], [-std]*2, 'g--', lw=3)

judge = det['RGB AGB Ratio'] >= 4
det = det[judge].reset_index(drop=True)
TRGB = np.array(det['TRGB'])
res = TRGB - np.mean(TRGB)
std = np.std(res)/c4(len(res))
ax[1].plot([3]*len(res), res, ls='', marker='D', ms=16, c='orange')
ax[1].plot([2.5, 3.5], [std]*2, 'g--', lw=3)
ax[1].plot([2.5, 3.5], [-std]*2, 'g--', lw=3)

det = pd.read_csv('../detection/ghosts_detection_v4.3.csv')
judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
judge = det['RGB AGB Ratio'] >= 4
det = det[judge].reset_index(drop=True)
TRGB = np.array(det['TRGB'])
res = TRGB - np.mean(TRGB)
std = np.std(res)/c4(len(res))
ax[1].plot([5]*len(res), res, ls='', marker='D', ms=16, c='orange')
ax[1].plot([4.5, 5.5], [std]*2, 'g--', lw=3)
ax[1].plot([4.5, 5.5], [-std]*2, 'g--', lw=3)

det = pd.read_csv('../detection/ghosts_detection_v4.4.csv')
judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
judge = det['RGB AGB Ratio'] >= 4
det = det[judge].reset_index(drop=True)
TRGB = np.array(det['TRGB'])
res = TRGB - np.mean(TRGB)
std = np.std(res)/c4(len(res))
ax[1].plot([7]*len(res), res, ls='', marker='D', ms=16, c='orange')
ax[1].plot([6.5, 7.5], [std]*2, 'g--', lw=3)
ax[1].plot([6.5, 7.5], [-std]*2, 'g--', lw=3)

judge = det['# star below tip'] >= 200
det = det[judge].reset_index(drop=True)
TRGB = np.array(det['TRGB'])
res = TRGB - np.mean(TRGB)
std = np.std(res)/c4(len(res))
ax[1].plot([9]*len(res), res, ls='', marker='D', ms=16, c='orange')
ax[1].plot([8.5, 9.5], [std]*2, 'g--', lw=3)
ax[1].plot([8.5, 9.5], [-std]*2, 'g--', lw=3)

det = pd.read_csv('../detection/ghosts_detection_v4.1.csv')
judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
judge = (det['RGB AGB Ratio'] >= 4) & (det['# star below tip'] >= 200)
det = det[judge].reset_index(drop=True)
TRGB = np.array(det['TRGB'])
res = TRGB - np.mean(TRGB)
std = np.std(res)/c4(len(res))
ax[1].plot([11]*len(res), res, ls='', marker='D', ms=16, c='orange')
ax[1].plot([10.5, 11.5], [std]*2, 'g--', lw=3)
ax[1].plot([10.5, 11.5], [-std]*2, 'g--', lw=3)

ax[1].set_xticks(np.arange(1, 12, 2))
ax[1].set_xticklabels(['raw', '$R$ cut', 'spatial', 'color band', '$N_{-,1.0}$ cut', 'smoothing'],fontsize=22)
ax[1].set_ylabel('Residual to Mean (mag)')
ax[1].set_ylim(-3.5, 3.5)
ax[1].set_yticks(np.arange(-3, 4, 1))
ax[1].grid()

#plt.savefig('Knobs.png')
plt.show()