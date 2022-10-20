#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 18:51:11 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MultipleLocator
from astropy.io import fits
from astropy.table import Table
from scipy import linalg
import itertools
import os

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

#change detection parameter here
info = pd.read_csv('ghosts_info_v3.csv')
ratio_min, nbt_min, peak_min = 3, 0, 0.6
sp_clip = True
smooth, weighting = 0.1, 'hatt'

agg_trgb = np.array([])
agg_ratio = np.array([])
agg_nbt = np.array([])
agg_field = np.array([])
agg_slope = np.array([])
agg_color = np.array([])
agg_raerr = np.array([])
agg_ext = np.array([])
agg_ratio2 = np.array([])

for i in range(len(info)):
    print(info['field'][i])
    #don't forget to change filename here
    data = pd.read_csv('clipped_csv/{:s}_10p.csv'.format(info['field'][i]))
    color, mag = data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS']
    band = np.array([info['slope_ac'][i], info['inter_ac'][i]])
    color, mag = selection(color, mag, band)
    mag, binval, edres = prod_edr(mag, smooth, weighting)
    trgb, ratio, nbt = findpeak(binval, edres, mag, ratio_min, nbt_min, peak_min)
    #change selection rules here, remember to change filename above at line 139
    judge = (ratio >= 4) & (nbt >= 100)
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
data.to_csv('ghosts_detection_v3.4.csv', index = False)