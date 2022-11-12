#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 18:37:39 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import optimize

def func(x, p0, p1):
    return p0 * x + p1

def error(p, x, y, err):
    return (func(p, x) - y) / err
    #return func(p, x) - y

plt.style.use('custom.mplstyle')
mpl.rcParams['text.usetex'] = False

#change version here
version, version_det = 4, 4.1

info = pd.read_csv('../info/ghosts_info_v{:s}.csv'.format(str(version)))
data = pd.read_csv('../detection/ghosts_detection_v{:s}.csv'.format(str(version_det)))
data['TRGB'] = data['TRGB'] - 10 * data['Extinction'] - data['A814']
judge = (data['RGB AGB Ratio'] >= 3) & (data['# star below tip'] >= 200) #change R,N limit here
data = data[judge].reset_index(drop=True)
gallist = ['0253','0891','2403','3031','4244','4565','4736','5236','7793','7814']

fig, ax = plt.subplots(figsize = (20, 8))
res, ratio, num = np.array([]), np.array([]), np.array([])

for i in range(len(gallist)):
    judge = [gallist[i] in data['field'][j] for j in range(len(data))]
    subdata = data[judge].reset_index(drop=True)
    TRGB, ra, n = [], [], []
    for j in range(len(subdata)):
        TRGB.append(subdata['TRGB'][j])
        ra.append(subdata['RGB AGB Ratio'][j])
        n.append(subdata['# star below tip'][j])
    if len(TRGB) > 1:
        TRGB = np.array(TRGB) - np.mean(np.array(TRGB))
        res = np.append(res, TRGB)
        ratio = np.append(ratio, ra)
        num  = np.append(num, n)

err = []
for i in range(len(ratio)):
    temp = np.exp(1.5*(3-ratio[i])) / (np.exp(1.5*(3-ratio[i]))+1) * (1/(num[i]-100))**0.1 * 2
    temp = temp ** 2 + 0.04 ** 2
    err.append(np.sqrt(temp))
err = np.array(err)
ax.errorbar(ratio, res, err, fmt='o', linewidth=2, markersize = 10, capsize=0)

#for removal of outlier
#judge = res>-0.2
#ratio, res, err = ratio[judge], res[judge], err[judge]

p0= [1, 20]
pfit, pcov = optimize.curve_fit(func, ratio, res, p0=p0, sigma=err, absolute_sigma=True)
perr = np.sqrt(np.diag(pcov))
print(perr)

k0, b = pfit[0], pfit[1]
unc = perr[0]
xx = np.array([2, 14])
yy = k0 * xx + b
ax.plot(xx, yy, color='orange', lw=3)

'''
# for calculation of chi^2/N
resp, ratiop, errp = res[res>-0.2], ratio[res>-0.2], err[res>-0.2]
residual = resp - (k0 * ratiop + b)
chi = np.zeros(len(residual))
for i in range(len(residual)):
    chi[i] = (residual[i] / errp[i]) ** 2
print(np.mean(chi))
'''

pfit1, pcov1 = optimize.curve_fit(func, ratio, res+err, p0=p0, sigma=err, absolute_sigma=True)
pfit2, pcov2 = optimize.curve_fit(func, ratio, res-err, p0=p0, sigma=err, absolute_sigma=True)
yy1 = (pfit1[0]) * xx + pfit1[1]
yy2 = (pfit2[0]) * xx + pfit2[1]
ax.fill_between(xx, yy1, yy2, color = 'orange', alpha=.3, linewidth=0)

ax.set_xlim(2, 12)
ax.set_xlabel('$R$')
ax.set_ylim(-0.5, 0.8)
ax.set_ylabel('Residual to Median (mag)')
ax.set_title('TRGB Tip-Contrast Relation')
ax.text(7.4, -0.4, 'slope={:.3f}$\pm${:.4f} mag/unit ratio'.format(k0, unc), color = 'red', fontsize=24)
ax.grid()

#plt.savefig('rtipscatter.png')
plt.show()
