#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 18:15:02 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import linalg

def selection(xx, yy, band):
    m, b = band[0], band[1]
    xx1 = (yy - b) / m
    xx2 = (yy - b) / m + 1.0
    judge = (xx > xx1) & (xx < xx2)
    xx, yy = xx[judge], yy[judge]
    return xx, yy

def gloess(mag, tau, bin_width = 0.01):
    fig_gloess, ax_gloess = plt.subplots(figsize = (10, 10))
    hist, bins, cont = ax_gloess.hist(mag, bins = np.arange(min(mag), max(mag) + bin_width, bin_width), color = 'black')
    plt.close()
    bin_centers = []
    for j in range(len(bins) - 1):
        bin_centers.append((bins[j] + bins[j + 1]) / 2)
    yest = np.zeros(len(hist))
    w = np.array([np.exp(- (bin_centers - bin_centers[i])**2/(2 * tau**2)) for i in range(len(hist))])
    for i in range(len(hist)):
        weights = w[i, :]
        b = np.array([np.sum(weights * hist), np.sum(weights * hist * bin_centers)])
        A = np.array([[np.sum(weights), np.sum(weights * bin_centers)],
                    [np.sum(weights * bin_centers), np.sum(weights * bin_centers * bin_centers)]])
        theta = linalg.solve(A, b)
        yest[i] = theta[0] + theta[1] * bin_centers[i]
    return yest, np.array(bin_centers[2: -2])

def poisson_sobel(hist):
    sobel = []
    for ii in range(2, len(hist) - 2):
        sobel.append(hist[ii - 1] * -1 * 1/np.sqrt(hist[ii - 1]) + hist[ii] * 0 + hist[ii + 1] * 1* 1/np.sqrt(hist[ii + 1]))
    sobel = np.array(sobel)
    return hist[2: -2], sobel

def hatt_sobel(hist):
    sobel = []
    for ii in range(2, len(hist) - 2):
        sobel.append((hist[ii + 1] - hist[ii - 1]) / np.sqrt(hist[ii + 1] + hist[ii - 1]))
    sobel = np.array(sobel)
    return hist[2: -2], sobel

def simple_sobel(hist):
    sobel = []
    for ii in range(2, len(hist) - 2):
        sobel.append(hist[ii-1] * -1 + hist[ii+1])
    sobel = np.array(sobel)
    return hist[2: -2], sobel

def prod_edr(mag, tau, weighting):
    hist, binval = gloess(mag, tau)
    if (weighting == 'poisson'):
        lumfunc, edres = poisson_sobel(hist)
    elif (weighting == 'hatt'):
        lumfunc, edres = hatt_sobel(hist)
    elif (weighting == 'simple'):
        lumfunc, edres = simple_sobel(hist)
    else:
        print('No such weighting type')
        return
    return mag, binval, edres

def sweep(yy):
    yy_min, yy_max = np.min(yy), np.max(yy)
    val = np.arange(np.round(yy_min * 10) / 10, np.round(yy_max * 10) / 10, 0.01)
    count = np.zeros(len(val))
    for i in range(len(val)):
        count[i] = sum((yy >= val[i] - 0.05) & (yy <= val[i] + 0.05))
    best = np.mean(val[count == np.max(count)])
    return best

color_list = ['#7e1e9c', '#15b01a', '#0343df', '#ff81c0', '#653700', '#e50000', '#95d0fc', '#029386', '#f97306', '#96f97b', '#c20078', '#ffff14',
              '#75bbfd', '#929591', '#89fe05', '#bf77f6', '#9a0eea', '#033500', '#06c2ac', '#c79fef', '#00035b', '#d1b26f', '#00ffff', '#13eac9',
              '#06470c', '#ae7181', '#35063e', '#01ff07', '#650021', '#6e750e', '#ff796c', '#e6daa6', '#0504aa', '#001146', '#cea2fd', '#000000',
              '#ff028d', '#ad8150', '#c7fdb5', '#ffb07c', '#677a04', '#cb416b', '#8e82fe', '#53fca1', '#aaff32', '#380282', '#ceb301', '#ffd1df',
              '#cf6275', '#0165fc', '#0cff0c', '#c04e01', '#04d8b2', '#01153e', '#3f9b0b', '#d0fefe', '#840000', '#be03fd', '#c0fb2d', '#a2cffe',
              '#dbb40c', '#8fff9f', '#580f41', '#4b006e', '#8f1402', '#014d4e', '#610023', '#aaa662', '#137e6d', '#7af9ab', '#02ab2e', '#9aae07']
name = ['01','02','03','04','05','06','07','08','09','10','11']

plt.style.use('custom.mplstyle')
mpl.rcParams['text.usetex'] = False

#change version here
version, version1, version2 = 4, 4.2, 4.1
info = pd.read_csv('../info/ghosts_info_v{:s}.csv'.format(str(version)))
judge = ['ngc3031' in info['field'][i] for i in range(len(info))]
info = info[judge].reset_index(drop=True)

fig = plt.figure(figsize = (42, 20))
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 3)

det = pd.read_csv('../detection/ghosts_detection_v{:s}.csv'.format(str(version1)))
judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
judge = (det['RGB AGB Ratio'] >= 3.5) & (det['# star below tip'] >= 50)
det = det[judge].reset_index(drop=True)
fit_trgb = sweep(np.array(det['TRGB']))

for i in range(len(info)):
    ax1.plot([i+1,i+1], [20, 35], 'k-', lw = 2)
ax1.plot([0,len(info)], [fit_trgb]*2, 'g--', lw=3, label = 'best fit')
ax1.legend(loc = 'lower center', ncol=2)

for i in range(len(info)):
    print(info['field'][i])
    filename = '../clipped_csv/{:s}_v{:s}.csv'.format(info['field'][i], str(version))
    data = pd.read_csv(filename)
    color, mag = data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS']
    band = [info['slope_ac'][i], info['inter_ac'][i]]
    
    color, mag = selection(color, mag, band)
    mag, binval, edres = prod_edr(mag, 0.05, 'hatt')
    ax2.plot(binval, edres, c=color_list[i], ls='-', lw=3, label='Halo Field{:s}'.format(name[i]))
    
    judge = [det['field'][j] == info['field'][i] for j in range(len(det))]
    subdet = det[judge].reset_index(drop=True)
    if len(subdet) == 0:
        continue
    
    ax1.plot([i+0.5]*len(subdet), subdet['TRGB'], c = '#faba6b', marker = 'D', markersize = 16, ls = '--', lw = 2)
    for j in range(len(subdet)):
        ax2.plot([subdet['TRGB'][j]]*2, [0, 10], c=color_list[i], ls='-', lw=2)
        if (subdet['RGB AGB Ratio'][j] >= 3.5) & (subdet['TRGB'][j] < 24.5) & (subdet['TRGB'][j] > 23.5):
            ax2.text(subdet['TRGB'][j], 0.3-0.02*i, '$R={:.1f}$'.format(subdet['RGB AGB Ratio'][j]), c = color_list[i], fontsize=24)
    
    judge = subdet['RGB AGB Ratio'] == np.max(np.array(subdet['RGB AGB Ratio']))
    subdet = subdet[judge].reset_index(drop=True)
    trgb = subdet['TRGB'][0]
    ax1.plot(i+0.5, trgb, color = '#fd7133', marker = 'D', markersize = 16, ls = '')

ax1.set_xlim(0, len(info))
ax1.set_ylim(25, 22.5)
ax1.set_xticks(np.arange(0.5, len(info), 1))
ax1.set_xticklabels(name)
ax1.set_xlabel('Halo Field')
ax1.set_ylabel('Tip Mag')
ax1.set_title('NGC 3031 smoothing = 0.05')

ax2.set_xlim(23.5, 24.3)
ax2.set_xlabel('F814W Mag')
ax2.set_ylim(0, 0.3)
ax2.set_ylabel('EDR')
ax2.legend(bbox_to_anchor=(1.1, -0.15), loc = 'upper center', ncol = 6)

ax1 = fig.add_subplot(2, 2, 2)
ax2 = fig.add_subplot(2, 2, 4)

det = pd.read_csv('../detection/ghosts_detection_v{:s}.csv'.format(str(version2)))
judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
judge = (det['RGB AGB Ratio'] >= 3.5) & (det['# star below tip'] >= 50)
det = det[judge].reset_index(drop=True)
fit_trgb = sweep(np.array(det['TRGB']))

for i in range(len(info)):
    ax1.plot([i+1,i+1], [20, 35], 'k-', lw = 2)
ax1.plot([0,len(info)], [fit_trgb]*2, 'g--', lw=3, label = 'best fit')
ax1.legend(loc = 'lower center', ncol=2)

for i in range(len(info)):
    print(info['field'][i])
    filename = '../clipped_csv/{:s}_v{:s}.csv'.format(info['field'][i], str(version))
    data = pd.read_csv(filename)
    color, mag = data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS']
    band = [info['slope_ac'][i], info['inter_ac'][i]]
    
    color, mag = selection(color, mag, band)
    mag, binval, edres = prod_edr(mag, 0.10, 'hatt')
    ax2.plot(binval, edres, c=color_list[i], ls='-', lw=3, label='Halo Field{:s}'.format(name[i]))
    
    judge = [det['field'][j] == info['field'][i] for j in range(len(det))]
    subdet = det[judge].reset_index(drop=True)
    if len(subdet) == 0:
        continue
    
    ax1.plot([i+0.5]*len(subdet), subdet['TRGB'], c = '#faba6b', marker = 'D', markersize = 16, ls = '--', lw = 2)
    for j in range(len(subdet)):
        ax2.plot([subdet['TRGB'][j]]*2, [0, 10], c=color_list[i], ls='-', lw=2)
        if (subdet['RGB AGB Ratio'][j] >= 3.5) & (subdet['TRGB'][j] < 24.5) & (subdet['TRGB'][j] > 23.5):
            ax2.text(subdet['TRGB'][j], 0.3-0.02*i, '$R={:.1f}$'.format(subdet['RGB AGB Ratio'][j]), c = color_list[i], fontsize=24)
    
    judge = subdet['RGB AGB Ratio'] == np.max(np.array(subdet['RGB AGB Ratio']))
    subdet = subdet[judge].reset_index(drop=True)
    trgb = subdet['TRGB'][0]
    ax1.plot(i+0.5, trgb, color = '#fd7133', marker = 'D', markersize = 16, ls = '')

ax1.set_xlim(0, len(info))
ax1.set_ylim(25, 22.5)
ax1.set_xticks(np.arange(0.5, len(info), 1))
ax1.set_xticklabels(name)
ax1.set_xlabel('Halo Field')
ax1.set_ylabel('Tip Mag')
ax1.set_title('NGC 3031 smoothing = 0.10')

ax2.set_xlim(23.5, 24.3)
ax2.set_xlabel('F814W Mag')
ax2.set_ylim(0, 0.3)
ax2.set_ylabel('EDR')

#plt.savefig('EDR.png')
plt.show()