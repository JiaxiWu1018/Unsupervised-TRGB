#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 18:41:11 2022

@author: michaelwu
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd

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

plt.style.use('custom.mplstyle')
mpl.rcParams['text.usetex'] = False

#change version here
version_ghosts, version_sn = 4, 1

info = pd.read_csv('../info/ghosts_info_v{:s}.csv'.format(str(version_ghosts)))

ra, pop, upp, col, slo, wid = [], [], [], [], [], []
for i in range(len(info)):
    data = pd.read_csv('../csv/{:s}.csv'.format(info['field'][i]))
    color, mag = np.array(data['MAG1_ACS'] - data['MAG2_ACS']), np.array(data['MAG2_ACS'])
    jcol, jmag = selection(color, mag, [info['slope_bc'][i], info['inter_bc'][i]])
    eddtip = info['eddtip'][i]
    
    RGB = np.sum((mag>=eddtip)&(mag<=eddtip+0.5))
    AGB = max(np.sum((mag<=eddtip)&(mag>=eddtip-0.5)),1)
    nbt = max(np.sum((mag>=eddtip)&(mag<=eddtip+1)),1)
    tipcolor = jcol[(jmag>=eddtip-0.1)&(jmag<=eddtip+0.1)]
    tipcolorm1 = jcol[(jmag>=eddtip+0.9)&(jmag<=eddtip+1.1)]
    tipcolorfull = color[(mag>=eddtip-0.1)&(mag<=eddtip+0.1)]
    
    judge_red, judge_blue = blue_selection(color, mag, info['red_lim'][i], info['faint_lim'][i])
    bcol, bmag = color[judge_blue], mag[judge_blue]
    bpop = np.sum((bmag>=eddtip-0.5)&(bmag<=eddtip+0.5))
    
    if (len(tipcolor) > 0) & (len(tipcolorm1) > 0):
        ra.append(RGB/AGB)
        pop.append(nbt)
        upp.append(bpop/nbt)
        col.append(np.median(tipcolor))
        slo.append(1/(np.median(tipcolorm1)-np.median(tipcolor)))
        wid.append(1.48*np.median(np.abs(tipcolorfull-np.median(tipcolor))))
        
ra, pop, upp, col, slo, wid = np.array(ra), np.array(pop), np.array(upp), np.array(col), np.array(slo), np.array(wid)
length = len(ra)

numra, binra, p1 = plt.hist(ra, bins = [0, 3, 4, 5, 7, 10, 1e2])
plt.close()
numra11 = numra / length
numpop, binpop, p1 = plt.hist(pop, bins = [0, 50, 100, 200, 500, 1000, 1e6])
plt.close()
numpop11 = numpop / length
numupp, binupp, p1 = plt.hist(upp, bins = [0, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 1])
plt.close()
numupp11 = numupp / length
numcol, bincol, p1 = plt.hist(col, bins = [0, 1.0, 1.2, 1.4, 1.6, 1.8, 5])
plt.close()
numcol11 = numcol / length
numslo, binslo, p1 = plt.hist(slo, bins = [-1e5, -8, -6, -4.5, -3, -2, 1e5])
plt.close()
numslo11 = numslo / length
numwid, binwid, p1 = plt.hist(wid, bins = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 10])
plt.close()
numwid11 = numwid / length

judge = (ra >= 3) & (pop >= 50)
ra, pop, upp, col, slo, wid = ra[judge], pop[judge], upp[judge], col[judge], slo[judge], wid[judge]

numra, binra, p1 = plt.hist(ra, bins = [0, 3, 4, 5, 7, 10, 1e2])
plt.close()
numra12 = numra / length
numpop, binpop, p1 = plt.hist(pop, bins = [0, 50, 100, 200, 500, 1000, 1e6])
plt.close()
numpop12 = numpop / length
numupp, binupp, p1 = plt.hist(upp, bins = [0, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 1])
plt.close()
numupp12 = numupp / length
numcol, bincol, p1 = plt.hist(col, bins = [0, 1.0, 1.2, 1.4, 1.6, 1.8, 5])
plt.close()
numcol12 = numcol / length
numslo, binslo, p1 = plt.hist(slo, bins = [-1e5, -8, -6, -4.5, -3, -2, 1e5])
plt.close()
numslo12 = numslo / length
numwid, binwid, p1 = plt.hist(wid, bins = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 10])
plt.close()
numwid12 = numwid / length

info = pd.read_csv('../../EDD/sn_info_v{:s}.csv'.format(str(version_sn)))

ra, pop, upp, col, slo, wid = [], [], [], [], [], []
for i in range(len(info)):
    data = pd.read_csv('../../EDD/csv/{:s}.csv'.format(info['field'][i]))
    color, mag = np.array(data['inst_vega_mag1'] - data['inst_vega_mag2']), np.array(data['inst_vega_mag2'])
    jcol, jmag = selection(color, mag, [info['slope_bc'][i], info['inter_bc'][i]])
    eddtip = info['eddtip'][i]
    
    RGB = np.sum((mag>=eddtip)&(mag<=eddtip+0.5))
    AGB = max(np.sum((mag<=eddtip)&(mag>=eddtip-0.5)),1)
    nbt = max(np.sum((mag>=eddtip)&(mag<=eddtip+1)),1)
    tipcolor = jcol[(jmag>=eddtip-0.1)&(jmag<=eddtip+0.1)]
    tipcolorm1 = jcol[(jmag>=eddtip+0.9)&(jmag<=eddtip+1.1)]
    tipcolorfull = color[(mag>=eddtip-0.1)&(mag<=eddtip+0.1)]
    
    judge_red, judge_blue = blue_selection(color, mag, info['red_lim'][i], info['faint_lim'][i])
    bcol, bmag = color[judge_blue], mag[judge_blue]
    bpop = np.sum((bmag>=eddtip-0.5)&(bmag<=eddtip+0.5))
    
    if (len(tipcolor) > 0) & (len(tipcolorm1) > 0):
        ra.append(RGB/AGB)
        pop.append(nbt)
        upp.append(bpop/nbt)
        col.append(np.median(tipcolor))
        slo.append(1/(np.median(tipcolorm1)-np.median(tipcolor)))
        wid.append(1.48*np.median(np.abs(tipcolorfull-np.median(tipcolor))))

ra, pop, upp, col, slo, wid = np.array(ra), np.array(pop), np.array(upp), np.array(col), np.array(slo), np.array(wid)
length = len(ra)

numra, binra, p1 = plt.hist(ra, bins = [0, 3, 4, 5, 7, 10, 1e2])
plt.close()
numra21 = numra / length
numpop, binpop, p1 = plt.hist(pop, bins = [0, 50, 100, 200, 500, 1000, 1e6])
plt.close()
numpop21 = numpop / length
numupp, binupp, p1 = plt.hist(upp, bins = [0, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 1])
plt.close()
numupp21 = numupp / length
numcol, bincol, p1 = plt.hist(col, bins = [0, 1.0, 1.2, 1.4, 1.6, 1.8, 5])
plt.close()
numcol21 = numcol / length
numslo, binslo, p1 = plt.hist(slo, bins = [-1e5, -8, -6, -4.5, -3, -2, 1e5])
plt.close()
numslo21 = numslo / length
numwid, binwid, p1 = plt.hist(wid, bins = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 10])
plt.close()
numwid21 = numwid / length

judge = (ra >= 3) & (pop >= 50)
ra, pop, upp, col, slo, wid = ra[judge], pop[judge], upp[judge], col[judge], slo[judge], wid[judge]
numra, binra, p1 = plt.hist(ra, bins = [0, 3, 4, 5, 7, 10, 1e2])
plt.close()
numra22 = numra / length
numpop, binpop, p1 = plt.hist(pop, bins = [0, 50, 100, 200, 500, 1000, 1e6])
plt.close()
numpop22 = numpop / length
numupp, binupp, p1 = plt.hist(upp, bins = [0, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 1])
plt.close()
numupp22 = numupp / length
numcol, bincol, p1 = plt.hist(col, bins = [0, 1.0, 1.2, 1.4, 1.6, 1.8, 5])
plt.close()
numcol22 = numcol / length
numslo, binslo, p1 = plt.hist(slo, bins = [-1e5, -8, -6, -4.5, -3, -2, 1e5])
plt.close()
numslo22 = numslo / length
numwid, binwid, p1 = plt.hist(wid, bins = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 10])
plt.close()
numwid22 = numwid / length

fig, ax = plt.subplots(1, 2, figsize=(36, 10))
ax[0].grid()
ax[1].grid()
xx = np.arange(6)
width = 0.35

ax[0].bar(xx-width/2, numra11, width, color='#1f77b4', alpha = 0.5, edgecolor='none', linestyle='--')
ax[1].bar(xx-width/2, numpop11, width, color='#1f77b4', alpha = 0.5, edgecolor='none', linestyle='--')

ax[0].bar(xx-width/2, numra12, width, color='#1f77b4', label='GHOSTS')
ax[1].bar(xx-width/2, numpop12, width, color='#1f77b4', label='GHOSTS')

ax[0].bar(xx+width/2, numra21, width, color='#ff7f0e', alpha = 0.5, edgecolor='none', linestyle='--')
ax[1].bar(xx+width/2, numpop21, width, color='#ff7f0e', alpha = 0.5, edgecolor='none', linestyle='--')

ax[0].bar(xx+width/2, numra22, width, color='#ff7f0e', label='SN Host')
ax[1].bar(xx+width/2, numpop22, width, color='#ff7f0e', label='SN Host')

label1 = ['$<3$', '$3 \sim 4$', '$4 \sim 5$', '$5 \sim 7$', '$7 \sim 10$', '$>10$']
label2 = ['$<50$', '$50 \sim 100$', '$100 \sim 200$', '$200 \sim 500$', '$500 \sim 1000$', '$>1000$']
label3 = ['$<10^{-4}$', '$10^{-4} \sim 5·10^{-4}$', '$5·10^{-4} \sim 10^{-3}$', '$10^{-3} \sim 5·10^{-3}$', '$5·10^{-3} \sim 10^{-2}$', '$>10^{-2}$']
label4 = ['$<1.0$', '$1.0 \sim 1.2$', '$1.2 \sim 1.4$', '$1.4 \sim 1.6$', '$1.6 \sim 1.8$', '$>1.8$']
label5 = ['$<-8$', '$-8 \sim -6$', '$-6 \sim -4.5$', '$-4.5 \sim -3$', '$-3 \sim -2$', '$>-2$']
label6 = ['$<0.1$', '$0.1 \sim 0.2$', '$0.2 \sim 0.4$', '$0.4 \sim 0.6$', '$0.6 \sim 0.8$', '$>0.8$']

ax[0].set_xticks(xx)
ax[0].set_xticklabels(label1)
ax[0].legend(loc='upper right')
ax[0].set_title('$R$')
ax[0].set_ylabel('Percentage of Fields')
ax[0].set_ylim(0, 0.5)
ax[1].set_xticks(xx)
ax[1].set_xticklabels(label2)
ax[1].legend(loc='upper left')
ax[1].set_title('$N_{-,1.0}$')
ax[1].set_ylim(0, 0.8)

#plt.savefig('galaxies.png')
plt.show()