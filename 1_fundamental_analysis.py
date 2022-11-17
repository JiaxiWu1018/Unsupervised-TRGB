#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 14:28:15 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
import os

def read_fits(filename):
    hdulist = fits.open(filename)
    dat = Table.read(hdulist[1], format = 'fits')
    photDF = dat.to_pandas()
    return photDF

flist = np.array(os.listdir('{:s}/fits'.format(os.getcwd())))
flist_new = []
for i in range(len(flist)):
    filename = flist[i]
    if 'fits' in filename:
        flist_new.append(filename)
flist = np.sort(np.array(flist_new))
gallist = np.array(['0253','0891','2403','3031','4244','4565','4736','5236','7793','7814'])
tip = np.array([23.82,26.04,23.50,23.93,24.13,26.39,24.24,24.52,23.76,26.81])

gal, field, eddtip, ra, dec = [], [], [], [], []
for i in range(len(flist)):
    filename = flist[i]
    gal.append(filename[:7])
    field.append(filename[:15])
    eddtip.append(tip[np.where(gallist == filename[3:7])[0][0]])
    
    data = read_fits('fits/{:s}'.format(filename))
    data.to_csv('csv/{:s}.csv'.format(field[-1]), index=False)
    
    ra.append(np.mean(data['RA']))
    dec.append(np.mean(data['DEC']))

radec = np.array(ra+dec).reshape(len(flist), 2)
np.savetxt('coordinate.txt', radec)
# GO TO https://irsa.ipac.caltech.edu/applications/DUST/ AND GENERATE ghosts extinction.txt TO CONTINUE

f = open('ghosts_extinction.txt', 'r')
result = []
colname = f.readline().strip().split()
for line in f.readlines():
    line = line.strip().split()
    for i in range(len(line)):
        line[i] = float(line[i])
    result.append(line)
f.close()
result = pd.DataFrame(result).dropna(axis = 0)
result.columns = colname
EBV = np.array(result['mean_E_B_V_SandF'])

info = pd.DataFrame()
info['galaxy'], info['field'], info['eddtip'] = np.array(gal), np.array(field), np.array(eddtip)
info['ra'], info['dec'] = np.array(ra), np.array(dec)
info['EBV'], info['A606'], info['A814'] = EBV, 2.285*EBV, 1.487*EBV

c_ra = np.array([11.888,35.63922,114.21417,148.88822,184.37358,189.08658,192.72109,204.25396,359.45763,0.81204])
c_dec = np.array([-25.28822,42.34915,65.60256,69.06529,37.80711,25.98767,41.12046,-29.86542,-32.59103,16.14542])
R = np.array([3.6,9.3,3.2,4.1,4.5,12.9,5.4,4.0,3.7,15.9])
posang = np.array([52,25,122,149,46,134,85,45,98,135])
inc1 = np.array([76,76,52,59,73,71,38,32,55,60])
t = (posang-90)*np.pi/180
inc = inc1*np.pi/180

gal = ['0253','0891','2403','3031','4244','4565','4736','5236','7793','7814']
dist, incl, posiang = [], [], []
for i in range(len(gal)):
    judge = [gal[i] in info['field'][j] for j in range(len(info))]
    flist = info[judge].reset_index(drop=True)
    ra, dec = flist['ra'], flist['dec']
    for j in range(len(flist)):
        deltara = -3600 * (ra[j]-c_ra[i]) * np.cos(dec[j]*np.pi/180)
        deltadec = (dec[j]-c_dec[i]) * 3600
        dis = (deltara*np.cos(t[i])+deltadec*np.sin(t[i])) ** 2 + (deltadec*np.cos(t[i])-deltara*np.sin(t[i]))**2/(np.cos(inc[i])**2)
        dis = np.sqrt(dis)
        diskpc = dis / 3600 * 3.14 / 180 * R[i] * 1000
        dist.append(diskpc)
        incl.append(inc1[i])
        posiang.append(posang[i])
dist, incl, posiang = np.array(dist), np.array(incl), np.array(posiang)
ext = 0.00414 * (dist / 70) ** (-0.84)

info['Inc'] = incl
info['Posang'] = posiang
info['Dist2nuc'] = dist
info['Ext'] = ext
info.to_csv('ghosts_analysis.csv', index=False)
