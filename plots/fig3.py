#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 17:34:13 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from scipy import linalg

def selection(xx, yy, band):
    m, b = band[0], band[1]
    xx1 = (yy - b) / m
    xx2 = (yy - b) / m + 1.0
    judge = (xx > xx1) & (xx < xx2)
    xx, yy = xx[judge], yy[judge]
    return xx, yy

def gloess_fig3(mag, tau, say=False, bin_width = 0.01):
    fig_gloess, ax_gloess = plt.subplots(figsize = (10, 10))
    hist, bins, cont = ax_gloess.hist(mag, bins = np.arange(min(mag), max(mag) + bin_width, bin_width), color = 'black')
    plt.close()
    bin_centers = []
    for j in range(len(bins) - 1):
        bin_centers.append((bins[j] + bins[j + 1]) / 2)
    
    if (say == True):
        ax2.plot(hist, bin_centers, c='grey', ls='-', lw=1)
    
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

def hatt_sobel(hist):
    sobel = []
    for ii in range(2, len(hist) - 2):
        sobel.append((hist[ii + 1] - hist[ii - 1]) / np.sqrt(hist[ii + 1] + hist[ii - 1]))
    sobel = np.array(sobel)
    return hist[2: -2], sobel

def qualify(trgb, mag, ratio_min, nbt_min):
    RGB_pop, AGB_pop, nbt = np.zeros(len(trgb)), np.zeros(len(trgb)), np.zeros(len(trgb))
    for i in range(len(trgb)):
        RGB_pop[i] = sum((mag < trgb[i] + 0.5) & (mag > trgb[i]))
        AGB_pop[i] = max(sum((mag < trgb[i]) & (mag > trgb[i] - 0.5)), 1)
        nbt[i] = sum((mag > trgb[i]) & (mag < trgb[i]+1))
    judge = ((RGB_pop / AGB_pop) >= ratio_min) & (nbt >= nbt_min)
    trgb = trgb[judge]
    return trgb, RGB_pop / AGB_pop, nbt

def findpeak(binval, func, mag, ratio_min, nbt_min, peak_min = 0.6, fwhm_min = 0.1):
    pos, peak, fwhm, ratio, nbt, peakq = [], [], [], [], [], []
    for i in range(1, len(func) - 1):
        if (func[i] > func[i - 1]) & (func[i] > func[i + 1]):
            left_edge, right_edge = i, i
            while (func[left_edge] > func[left_edge - 1]) & (func[left_edge] > func[i] / 2) & (left_edge > 0):
                left_edge -= 1
            while (func[right_edge] > func[right_edge + 1]) & (func[right_edge] > func[i] / 2) & (right_edge < len(func) - 2):
                right_edge += 1
            
            test, RA, NBT = qualify(np.array([binval[i]]), mag, ratio_min, nbt_min)
            pos.append(i)
            ratio.append(RA[0])
            nbt.append(NBT[0])
            peak.append(func[i])
            fwhm.append(right_edge - left_edge)
            if len(test) == 1:
                peakq.append(func[i])
            i = right_edge
    pos_good, ratio_good, nbt_good = [], [], []
    if len(peakq) > 0:
        maxv = max(peakq)
        for i in range(len(pos)):
            if (peak[i] >= peak_min * maxv) & (fwhm[i] > fwhm_min):
                pos_good.append(pos[i])
                ratio_good.append(ratio[i])
                nbt_good.append(nbt[i])
    trgb = binval[pos_good]
    return trgb, np.array(ratio_good), np.array(nbt_good)

def blue_selection(xx, yy, xmin, ymax):
    judge_blue = (xx <= xmin) & (yy <= ymax)
    judge_red = [not i for i in judge_blue]
    return judge_red, judge_blue

def num_den(xx, yy, bins = 50):
    fig, ax = plt.subplots(figsize = (6, 6))
    hist, binx, biny, imag = ax.hist2d(xx, yy, bins = bins, cmap = cm.binary)
    plt.close()
    return hist, binx, biny

def gloess2d(xx, yy, zz, tau = 0.2):
    yest = np.zeros(zz.shape)
    xlen, ylen = np.arange(xx.shape[0])/10, np.arange(xx.shape[1])/10
    ylen, xlen = np.meshgrid(ylen, xlen)
    w = np.array([])
    for i in range(xx.shape[0]):
        for j in range(xx.shape[1]):
            w_new = np.array(np.exp(-((xlen - i/10) ** 2 + (ylen - j/10) ** 2) / (2 * tau ** 2)))
            w = np.append(w, w_new)
    w = w.reshape(xx.shape[0], xx.shape[1], xx.shape[0], xx.shape[1])
    for i in range(zz.shape[0]):
        for j in range(zz.shape[1]):
            weight = w[i, j, :, :]
            A = np.array([[np.sum(weight*xx**2),np.sum(weight*xx*yy),np.sum(weight*xx)],
                          [np.sum(weight*xx*yy),np.sum(weight*yy**2),np.sum(weight*yy)],
                          [np.sum(weight*xx),np.sum(weight*yy),np.sum(weight)]])
            b = np.array([np.sum(weight*xx*zz), np.sum(weight*yy*zz), np.sum(weight*zz)])
            theta = linalg.solve(A, b)
            yest[i, j] = theta[0]*xx[i, j]+theta[1]*yy[i, j]+theta[2]
    return yest

def spatial_clip(den, binx, biny, data, spclip):
    if (spclip == 'no'):
        return data
    if (spclip[-1] == 'p'):
        trun = int(spclip[:-1]) / 100 * np.max(den)
    else:
        trun = float(spclip)
    for i in range(den.shape[0]):
        for j in range(den.shape[1]):
            if den[i, j] >= trun:
                judge = (data['X'] >= binx[i]) & (data['X'] <= binx[i + 1]) & (data['Y'] >= biny[j]) & (data['Y'] <= biny[j + 1])
                judge = [not item for item in judge]
                data = data[judge]
    return data

plt.style.use('custom.mplstyle')
mpl.rcParams['text.usetex'] = False

#change version and field here
version = 4
field = 'ngc3031-field05'
title = 'NGC3031 Halo Field04'

info = pd.read_csv('../info/ghosts_info_v{:s}.csv'.format(str(version)))
judge = [info['field'][i] == field for i in range(len(info))]
info = info[judge].reset_index(drop=True)
data = pd.read_csv('../csv/{:s}.csv'.format(field))
color, mag = data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS']

fig = plt.figure(figsize = (20, 11))
ax1 = plt.subplot(1, 10, (1, 4))
ax2 = plt.subplot(1, 10, (5, 7))
ax3 = plt.subplot(1, 10, (8, 10))
fig.suptitle(title, fontweight='bold')

m, b = info['slope_bc'][0], info['inter_bc'][0]
yy = np.array([np.ceil(info['faint_lim'][0]+1), np.floor(info['faint_lim'][0]-5)])
xlist = np.arange(0, 1.01, 0.01)
for i in range(len(xlist)):
    xx = (yy - b) / m + xlist[i]
    ax1.plot(xx, yy, c = '#6f94cd', alpha=0.3)

ax1.set_xlim(-1, 3.5)
ax1.set_ylim(info['faint_lim'][0]+1, info['faint_lim'][0]-5)
ax1.set_xlabel('Color: F606W-F814W')
ax1.set_ylabel('Brightness: F814W')
ax1.set_yticks(np.arange(np.ceil(info['faint_lim'][0]-5), np.ceil(info['faint_lim'][0]+1)))
ax1.set_title('Color Magnitude Diagram')

hist, binval = gloess_fig3(mag, tau=0.10)
lumfunc, edres = hatt_sobel(hist)
ax3.plot(edres, binval, c='orange', ls='-', lw=4, label = 'raw')

trgb, ratio, nbt = findpeak(binval, edres, mag, 3, 0, 0.6)

judge_red, judge_blue = blue_selection(data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS'], info['red_lim'][0], info['faint_lim'][0])
data_red, data_blue = data[judge_red].reset_index(drop=True), data[judge_blue].reset_index(drop=True)
Z_blue, binx, biny = num_den(np.array(data_blue['X']), np.array(data_blue['Y']))
Z_red, binx, biny = num_den(np.array(data_red['X']), np.array(data_red['Y']), (binx, biny))
X, Y = (binx[1:]+binx[:-1])/2, (biny[1:]+biny[:-1])/2
Y, X = np.meshgrid(Y, X)

Z_ratio = np.zeros(np.shape(Z_red))
for j in range(len(Z_blue)):
    for k in range(len(Z_blue[0, :])):
        if Z_blue[j, k] == 0:
            Z_ratio[j, k] = np.max(Z_red)
        else:
            Z_ratio[j, k] = Z_red[j, k] / Z_blue[j, k]

Z_smooth = gloess2d(X, Y, Z_blue)
data_clipped = spatial_clip(Z_smooth, binx, biny, data, '10p').reset_index(drop=False)
judge = [True] * len(data)
for i in range(len(data_clipped)):
    judge[data_clipped['index'][i]] = False
data_clip = data[judge].reset_index(drop=True)
color, mag = data_clip['MAG1_ACS']-data_clip['MAG2_ACS'], data_clip['MAG2_ACS']
ax1.plot(color, mag, c = 'red', marker = 'x', ls = '', markersize=3)
color, mag = data_clipped['MAG1_ACS']-data_clipped['MAG2_ACS'], data_clipped['MAG2_ACS']
ax1.plot(color, mag, 'k.', markersize=3)

color, mag = selection(color, mag, [m, b])
hist, binval = gloess_fig3(mag, tau=0.10, say=True)
lumfunc, edres = hatt_sobel(hist)
ax2.plot(lumfunc, binval, 'k-', lw = 4)
ax3.plot(edres, binval, 'k-', lw=4, label = 'clipped')

trgb, ratio, nbt = findpeak(binval, edres, mag, 3, 0, 0.6)
judge = (ratio >= 4)
trgb, ratio, nbt = trgb[judge], ratio[judge], nbt[judge]
for i in range(len(trgb)):
    ax1.plot([-1, 3.5], [trgb[i]]*2, 'r-', lw=3)
    ax2.plot([0, 1e3], [trgb[i]]*2, 'r-', lw=3)
    ax3.plot([0, 1], [trgb[i]]*2, 'r-', lw=3)

ax2.set_xlim(0, 45)
ax2.set_xlabel('# per bin')
ax2.set_ylim(info['faint_lim'][0]+1, info['faint_lim'][0]-5)
ax2.set_title('Luminosity Function')
ax2.set_yticks([])

ax3.legend(loc='upper left')
ax3.set_xlim(0.18, 0)
ax3.set_ylim(info['faint_lim'][0]+1, info['faint_lim'][0]-5)
ax3.set_yticks([])
ax3.set_title('Edge Detection Response')

#plt.savefig('explainer.png')
plt.show()