#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 19:00:47 2022

@author: michaelwu
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.style.use('custom.mplstyle')
mpl.rcParams['text.usetex'] = False

age = np.array([1e9,2e9,3e9,4e9,5e9,6e9,7e9,8e9,9e9,1e10])/1e9
ratio1 = [0.931707,4.65217,6.73979,7.92388,7.85154,8.37684,14.2609,16.8758,20.6722,20.1595]
ratio2 = [1.02500,7.30785,12.1772,16.9628,18.2526,18.9903,16.9549,15.7017,14.7965,15.1304]

fig, ax = plt.subplots(figsize=(16, 8))
ax.plot(age, ratio1, ls='-', c = 'orange', marker='s', lw=4, ms=12, label='[M/H]$=-2$')
ax.plot(age, ratio2, ls='-', marker='D', lw=4, ms=12, label='[M/H]$=-1.5$')
ax.set_xlabel('stellar pop age (Gyr)')
ax.set_ylabel('$R$')
ax.legend(loc='best', ncol=2)
ax.grid()
ax.set_title('Tip Ratio vs Age')

#plt.savefig('age.png')
plt.show()