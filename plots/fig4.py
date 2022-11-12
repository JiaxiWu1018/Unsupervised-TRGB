#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 17:46:37 2022

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

def color_band(color, mag, slope_range):
    num = 0
    for b in np.arange(25, 35, 0.1):
        for m in np.arange(slope_range[0], slope_range[1], 0.1):
            xx1 = (mag - b) / m
            xx2 = (mag - b) / m + 1.0
            judge = (color < xx2) & (color > xx1)
            if np.sum(judge) > num:
                num = np.sum(judge)
                slope, inter = m, b
    return np.array([slope, inter])

plt.style.use('custom.mplstyle')
mpl.rcParams['text.usetex'] = False

#change version here
version = 4
field = 'ngc3031-field10'
title = 'NGC3031 Halo Field 09'

vinfo = pd.read_csv('../version_code2.csv')
vinfo = vinfo[vinfo['version']==version].reset_index(drop=True)
spclip, slope_range = vinfo['spatial_clip'][0], [vinfo['slope_min'][0], vinfo['slope_max'][0]+0.05]

info = pd.read_csv('../info/ghosts_info_v{:s}.csv'.format(str(version)))
judge = [info['field'][i] == field for i in range(len(info))]
info = info[judge].reset_index(drop=True)

for i in range(len(info)):
    print(info['field'][i])
    data = pd.read_csv('../csv/{:s}.csv'.format(info['field'][i]))
    color, mag = data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS']
    blim = np.round(info['faint_lim'][i]*2)/2 + 1
    
    judge_red, judge_blue = blue_selection(data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS'], info['red_lim'][i], info['faint_lim'][i])
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
    data_clip = spatial_clip(Z_smooth, binx, biny, data, spclip).reset_index(drop=True)
    
    fig = plt.figure(figsize = (32, 20))
    ax1 = fig.add_subplot(2, 17, (1,5))
    cs = ax1.contourf(X, Y, Z_ratio, cmap = 'Reds')
    ax1.set_title('Red to Blue Ratio')
    cbar = fig.colorbar(cs)
    ax2 = fig.add_subplot(2, 17, (7,11))
    cs = ax2.contourf(X, Y, Z_blue, cmap = 'Blues')
    ax2.set_title('Blue Star Pop')
    cbar = fig.colorbar(cs)
    ax3 = fig.add_subplot(2, 17, (13,17))
    cs = ax3.contourf(X, Y, Z_smooth, cmap = 'Blues')
    ax3.set_title('Blue Star Pop Smoothed')
    cbar = fig.colorbar(cs)
    
    ax4 = fig.add_subplot(2, 18, (19, 22))
    ax4.plot(data_red['MAG1_ACS']-data_red['MAG2_ACS'], data_red['MAG2_ACS'], 'k.', markersize = 2)
    ax4.plot(data_blue['MAG1_ACS']-data_blue['MAG2_ACS'], data_blue['MAG2_ACS'], 'b.', markersize = 2)
    ax4.set_ylim(blim, blim-6)
    ax4.set_xlim(-1, 3.5)
    ax4.set_xlabel('F606W-F814W')
    ax4.set_ylabel('F814W')
    ax4.set_title('CMD raw')
    
    color, mag = data_clip['MAG1_ACS']-data_clip['MAG2_ACS'], data_clip['MAG2_ACS']
    judge_red, judge_blue = blue_selection(color, mag, info['red_lim'][i], info['faint_lim'][i])
    data_red, data_blue = data_clip[judge_red], data_clip[judge_blue]
    ax5 = fig.add_subplot(2, 18, (25, 30))
    ax5.plot(data_red['X'], data_red['Y'], 'r.', markersize = 2)
    ax5.plot(data_blue['X'], data_blue['Y'], 'b.', markersize = 2)
    ax5.set_xlabel('X')
    ax5.set_ylabel('Y')
    ax5.set_title('Spatial Distribution')
    
    ax6 = fig.add_subplot(2, 18, (33, 36))
    ax6.plot(data_clip['MAG1_ACS']-data_clip['MAG2_ACS'], data_clip['MAG2_ACS'], 'k.', markersize = 2)
    ax6.set_xlim(-1, 3.5)
    ax6.set_xlabel('F606W-F814W')
    ax6.set_ylim(blim, blim-6)
    ax6.set_ylabel('F814W')
    ax6.set_title('CMD clipped')
    
    ax2.text(-200, 4750, title, fontsize=36)
    ax5.text(200, 2200, 'In direction', fontsize=28)
    ax5.text(350, 1900, 'of disk', fontsize=28)
    
    #plt.savefig('spatial.png')
    plt.show()