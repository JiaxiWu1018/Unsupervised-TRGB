#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 15:13:28 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import linalg

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

#change spatial clipping trun here
def spatial_clip(den, binx, biny, data):
    trun = 0.1 * np.max(den)
    for i in range(den.shape[0]):
        for j in range(den.shape[1]):
            if den[i, j] >= trun:
                judge = (data['X'] >= binx[i]) & (data['X'] <= binx[i + 1]) & (data['Y'] >= biny[j]) & (data['Y'] <= biny[j + 1])
                judge = [not item for item in judge]
                data = data[judge]
    return data

#change searching scope here
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

info = pd.read_csv('ghosts_analysis.csv')
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

slope, inter = np.zeros(len(info)), np.zeros(len(info))
for i in range(len(info)):
    print(info['field'][i])
    data = pd.read_csv('csv/{:s}.csv'.format(info['field'][i]))
    color, mag = data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS']
    blim = np.round(np.max(mag)*2)/2
    red_lim = 0.3 + info['A606'][i] - info['A814'][i]
    faint_lim = info['inter_bc'][i] + info['slope_bc'][i] * red_lim
    
    judge_red, judge_blue = blue_selection(data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS'], red_lim, faint_lim)
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
    data_clip = spatial_clip(Z_smooth, binx, biny, data).reset_index(drop=True)
    
    fig = plt.figure(figsize = (16, 10))
    ax1 = fig.add_subplot(2, 3, 1)
    cs = ax1.contourf(X, Y, Z_ratio, cmap = 'Reds')
    ax1.set_title('Red to Blue Ratio')
    cbar = fig.colorbar(cs)
    ax2 = fig.add_subplot(2, 3, 2)
    cs = ax2.contourf(X, Y, Z_blue, cmap = 'Blues')
    ax2.set_title('Blue Star Pop')
    cbar = fig.colorbar(cs)
    ax3 = fig.add_subplot(2, 3, 3)
    cs = ax3.contourf(X, Y, Z_smooth, cmap = 'Blues')
    ax3.set_title('Blue Star Pop Smoothed')
    cbar = fig.colorbar(cs)
    
    ax4 = fig.add_subplot(2, 18, (20, 23))
    ax4.plot(data_red['MAG1_ACS']-data_red['MAG2_ACS'], data_red['MAG2_ACS'], 'k.', markersize = 1)
    ax4.plot(data_blue['MAG1_ACS']-data_blue['MAG2_ACS'], data_blue['MAG2_ACS'], 'b.', markersize = 1)
    #ax4.plot([red_lim]*2, [base, base-6], 'r-', lw=1)
    ax4.set_ylim(blim, blim-6)
    ax4.set_xlim(-1, 3.5)
    ax4.set_xlabel('F606-F814')
    ax4.set_ylabel('F814')
    ax4.set_title('CMD raw')
    
    color, mag = data_clip['MAG1_ACS']-data_clip['MAG2_ACS'], data_clip['MAG2_ACS']
    judge_red, judge_blue = blue_selection(color, mag, red_lim, faint_lim)
    data_red, data_blue = data_clip[judge_red], data_clip[judge_blue]
    ax5 = fig.add_subplot(2, 18, (25, 30))
    ax5.plot(data_red['X'], data_red['Y'], 'r.', markersize = 1)
    ax5.plot(data_blue['X'], data_blue['Y'], 'b.', markersize = 1)
    ax5.set_xlabel('X')
    ax5.set_ylabel('Y')
    ax5.set_title('Spatial Distribution')
    
    ax6 = fig.add_subplot(2, 18, (32, 35))
    ax6.plot(data_clip['MAG1_ACS']-data_clip['MAG2_ACS'], data_clip['MAG2_ACS'], 'k.', markersize = 1)
    ax6.set_xlim(-1, 3.5)
    ax6.set_xlabel('F606-F814')
    ax6.set_ylim(blim, blim-6)
    ax6.set_ylabel('F814')
    ax6.set_title('CMD clipped')
    plt.savefig('clipping_map/{:s}_10p.png'.format(info['field'][i]))
    plt.show()
    
    data_clip.to_csv('clipped_csv/{:s}_10p.csv'.format(info['field'][i]))
    slope[i], inter[i] = color_band(color, mag)
    
info['slope_ac'], info['inter_ac'] = slope, inter
info.to_csv('ghosts_info_v1.csv', index=False)

i = 0
while i < len(info):
    gal = info['galaxy'][i]
    judge = [(gal == info['galaxy'][j]) for j in range(len(info))]
    subinfo = info[judge].reset_index(drop=True)
    i += len(subinfo)
    
    fig, ax = plt.subplots(1, len(subinfo), figsize = (3.5*len(subinfo), 6))
    for j in range(len(subinfo)):
        data = pd.read_csv('clipped_csv/{:s}_10p.csv'.format(subinfo['field'][j]))
        color, mag = data['MAG1_ACS']-data['MAG2_ACS'], data['MAG2_ACS']
        ax[j].plot(color, mag, 'k.', ms=1)
        
        m, b = subinfo['slope_ac'][j], subinfo['inter_ac'][j]
        yy = np.array([35, 20])
        xx1, xx2 = (yy - b) / m, (yy - b) / m + 1.0
        ax[j].plot(xx1, yy, c='orange', ls='--', lw=2)
        ax[j].plot(xx2, yy, c='orange', ls='--', lw=2)
        
        ax[j].set_xlim(-1, 3.5)
        blim = np.round(np.max(mag)*2)/2
        ax[j].set_ylim(blim, blim-6)
        ax[j].set_title(subinfo['field'][j])
    plt.savefig('clipped_cmd/{:s}_10p.png'.format(gal))
    plt.show()
