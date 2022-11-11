#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 19:02:30 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import linalg

def c4(n):
    k = n//2
    if n%2 == 0:
        return np.sqrt(2/(np.pi*(2*k-1)))*(2**(2*k-2))*(math.factorial(k-1)**2)/math.factorial(2*k-2)
    else:
        return np.sqrt(np.pi/k)*math.factorial(2*k-1)/(2**(2*k-1))/(math.factorial(k-1)**2)

def color_band(color, mag, m, w):
    num = 0
    for b in np.arange(25, 35, 0.1):
        xx1 = (mag - b) / m
        xx2 = (mag - b) / m + w
        judge = (color < xx2) & (color > xx1)
        if np.sum(judge) > num:
            num = np.sum(judge)
            slope, inter = m, b
    return np.array([slope, inter])

def selection(xx, yy, band, w):
    m, b = band[0], band[1]
    xx1 = (yy - b) / m
    xx2 = (yy - b) / m + w
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

plt.style.use('custom.mplstyle')
mpl.rcParams['text.usetex'] = False

#change version and knobs for optimization here
version = 4
knob = 'Weighting' #Slope, Width, Smoothing, Weighting (capitalized)

info = pd.read_csv('../info/ghosts_info_v{:s}.csv'.format(str(version)))
ratio_min, nbt_min, peak_min = 3, 0, 0.6
slope, width, smooth, weighting = -6, 1, 0.10, 'hatt'
stdra, per, tpf = [], [], []

if (knob == 'Slope'):
    l = np.linspace(-10, -2, 17)
if (knob == 'Width'):
    l = np.linspace(0.5, 2, 16)
if (knob == 'Smoothing'):
    l = np.linspace(0.03, 0.12, 19)
if (knob == 'Weighting'):
    l = np.array(['simple', 'poisson', 'hatt'])

for j in range(len(l)):
    if (knob == 'Slope'):
        slope = l[j]
    if (knob == 'Width'):
        width = l[j]
    if (knob == 'Smoothing'):
        smooth = l[j]
    if (knob == 'Weighting'):
        weighting = l[j]
    print(l[j])
    
    agg_trgb = np.array([])
    agg_ratio = np.array([])
    agg_nbt = np.array([])
    agg_field = np.array([])
    agg_slope = np.array([])
    agg_color = np.array([])
    agg_raerr = np.array([])
    agg_ext = np.array([])
    agg_ratio2 = np.array([])
    agg_a814 = np.array([])
    agg_a606 = np.array([])

    for i in range(len(info)):
        data = pd.read_csv('../clipped_csv/{:s}_v4.csv'.format(info['field'][i]))
        color, mag = data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS']
        
        band = color_band(color, mag, slope, width)
        color, mag = selection(color, mag, band, width)
        mag, binval, edres = prod_edr(mag, smooth, weighting)
        trgb, ratio, nbt = findpeak(binval, edres, mag, ratio_min, nbt_min, peak_min)
        judge = (ratio >= 1.5) & (nbt >= 50)
        trgb, ratio, nbt = trgb[judge], ratio[judge], nbt[judge]
        
        tipcol, tipsl, raerr, ra2 = np.array([]), np.array([]), np.array([]), np.array([])
        judge = [True]*len(trgb)
        for k in range(len(trgb)):
            judgetip = (mag < trgb[k]+0.1) & (mag > trgb[k]-0.1)
            judge1b = (mag < trgb[k]+1.1) & (mag > trgb[k]+0.9)
            if (np.sum(judgetip) == 0):
                judge[k]=False
                continue
            elif (np.sum(judge1b) == 0):
                maxi = np.max(mag)
                judge1b = (mag > maxi-0.1)
            coltip = np.mean(color[judgetip])
            col1b = np.mean(color[judge1b])
            RGB = np.sum((mag>trgb[k])&(mag<trgb[k]+0.5))
            AGB = max(np.sum((mag<trgb[k])&(mag>trgb[k]-0.5)),1)
            tipra = RGB/AGB
            tipraerr = tipra * np.sqrt(1/RGB + 1/AGB)
            tipcol = np.append(tipcol, np.array([coltip]))
            tipsl = np.append(tipsl, np.array([1 / (coltip - col1b)]))
            raerr = np.append(raerr, np.array([tipraerr]))
            RGB2 = np.sum((mag>trgb[k])&(mag<trgb[k]+1))
            AGB2 = max(np.sum((mag<trgb[k])&(mag>trgb[k]-1)),1)
            ra2 = np.append(ra2, np.array([RGB2/AGB2]))
        trgb, ratio, nbt = trgb[judge], ratio[judge], nbt[judge]
        
        agg_trgb = np.append(agg_trgb, trgb)
        agg_ratio = np.append(agg_ratio, ratio)
        agg_nbt = np.append(agg_nbt, nbt)
        agg_field = np.append(agg_field, np.array([info['field'][i]]*len(trgb)))
        agg_ext = np.append(agg_ext, np.array([info['Ext'][i]]*len(trgb)))
        agg_a606 = np.append(agg_a606, np.array([info['A606'][i]]*len(trgb)))
        agg_a814 = np.append(agg_a814, np.array([info['A814'][i]]*len(trgb)))
        agg_slope = np.append(agg_slope, tipsl)
        agg_color = np.append(agg_color, tipcol)
        agg_raerr = np.append(agg_raerr, raerr)
        agg_ratio2 = np.append(agg_ratio2, ra2)

    agg_galaxy = np.array([agg_field[i][:7] for i in range(len(agg_field))])
    data = pd.DataFrame()
    data['galaxy'] = agg_galaxy
    data['field'] = agg_field
    data['TRGB'] = agg_trgb
    data['RGB AGB Ratio'] = agg_ratio
    data['Ratio Error'] = agg_raerr
    data['RGB AGB Ratio pm 1'] = agg_ratio2
    data['# star below tip'] = agg_nbt
    data['Slope'] = agg_slope
    data['Tip Color'] = agg_color
    data['Extinction'] = agg_ext
    data['A606'] = agg_a606
    data['A814'] = agg_a814
    judge = (data['RGB AGB Ratio'] >= 4) & (data['# star below tip'] >= 200)
    data = data[judge].reset_index(drop=True)
    
    gallist = ['0253','0891','2403','3031','4244','4565','4736','5236','7793','7814']
    std, num, field = [], [], []
    for k in range(len(data)):
        if data['field'][k] not in field:
            field.append(data['field'][k])
    for k in range(len(gallist)):
        gal = gallist[k]
        judge = [gal in data['field'][m] for m in range(len(data))]
        if np.sum(judge) >= 2:
            std.append(np.std(data['TRGB'][judge])/c4(int(np.sum(judge))))
            num.append(int(np.sum(judge)))

    summ = 0
    for k in range(len(std)):
        summ += std[k] * num[k]
    stdra.append(summ / np.sum(np.array(num)))
    per.append(100 * len(field) / 51)
    tpf.append(len(data)/len(field))

fig, ax = plt.subplots(3, 1, figsize=(18, 18))
l1 = False
if l[0] == 'simple':
    l = [1, 2, 3]
    l1 = True

ax[0].plot(l, stdra, marker='*', ms=15, lw=4)
ax[0].set_xlim(np.min(l)-0.1, np.max(l)+0.1)
ax[0].set_ylim(0.01, 0.5)
ax[0].set_yscale('log')
ax[0].set_yticks([1e-2,5e-2,1e-1,5e-1])
ax[0].set_yticklabels(['0.01','0.05', '0.10','0.50'])
ax[0].plot([np.min(l)-1, np.max(l)+1], [0.1]*2, c='grey', ls='--', lw=2)
ax[0].plot([np.min(l)-1, np.max(l)+1], [0.05]*2, c='grey', ls='--', lw=2)
ax[0].set_title('Optimization of Multifield TRGB Dispersion on {:s}'.format(knob))
ax[0].set_ylabel('$\sigma_{tot}$ (mag)')

ax[1].plot(l, per, marker='*', ms=15, lw=4)
ax[1].set_xlim(np.min(l)-0.1, np.max(l)+0.1)
ax[1].set_ylim(45, 85)
ax[1].plot([np.min(l)-1, np.max(l)+1], [70]*2, c='grey', ls='--', lw=2)
ax[1].plot([np.min(l)-1, np.max(l)+1], [50]*2, c='grey', ls='--', lw=2)
ax[1].set_ylabel('$P_{valid}$')

ax[2].plot(l, tpf, marker='*', ms=15, lw=4)
ax[2].set_xlim(np.min(l)-0.1, np.max(l)+0.1)
ax[2].set_ylim(0.95, 1.4)
ax[2].plot([np.min(l)-1, np.max(l)+1], [1.3]*2, c='grey', ls='--', lw=2)
ax[2].plot([np.min(l)-1, np.max(l)+1], [1.0]*2, c='grey', ls='--', lw=2)
ax[2].set_ylabel('$N_{TPF}$')
ax[2].set_xlabel(knob)
if l1 == True:
    ax[0].set_xticks([1, 2, 3])
    ax[0].set_xticklabels(['simple', 'poisson', 'hatt'])
    ax[1].set_xticks([1, 2, 3])
    ax[1].set_xticklabels(['simple', 'poisson', 'hatt'])
    ax[2].set_xticks([1, 2, 3])
    ax[2].set_xticklabels(['simple', 'poisson', 'hatt'])

#plt.savefig('optimization_{:s}'.format(knob))
plt.show()
