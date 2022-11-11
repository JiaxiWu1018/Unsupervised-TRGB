#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 18:31:41 2022

@author: michaelwu
"""

import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

def count_field(det):
    field, gal = [], []
    for i in range(len(det)):
        if det['field'][i] not in field:
            field.append(det['field'][i])
        if det['galaxy'][i] not in gal:
            gal.append(det['galaxy'][i])
    return field, gal

def c4(n):
    k = n//2
    if n%2 == 0:
        return np.sqrt(2/(np.pi*(2*k-1)))*(2**(2*k-2))*(math.factorial(k-1)**2)/math.factorial(2*k-2)
    else:
        return np.sqrt(np.pi/k)*math.factorial(2*k-1)/(2**(2*k-1))/(math.factorial(k-1)**2)

def cal_dispersion(gal, det):
    summ, length = 0, 0
    for i in range(len(gal)):
        judge = [det['galaxy'][j] == gal[i] for j in range(len(det))]
        subdet = det[judge].reset_index(drop=True)
        if len(subdet) >= 2:
            std = np.std(np.array(subdet['TRGB'])/c4(len(subdet)))
            summ += len(subdet) * std
            length += len(subdet)
    std = summ / length
    return std

det = pd.read_csv('../detection/ghosts_detection_v1.2.csv')
field, gal = count_field(det)
std = cal_dispersion(gal, det)
print('all', std, len(field)/50 * 100, len(det)/len(field))

judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
field, gal = count_field(det)
std = cal_dispersion(gal, det)
print('3031', std, len(field)/11 * 100, len(det)/len(field))

det = pd.read_csv('../detection/ghosts_detection_v4.3.csv')
field, gal = count_field(det)
std = cal_dispersion(gal, det)
print('all', std, len(field)/50 * 100, len(det)/len(field))

judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
field, gal = count_field(det)
std = cal_dispersion(gal, det)
print('3031', std, len(field)/11 * 100, len(det)/len(field))

det = pd.read_csv('../detection/ghosts_detection_v4.4.csv')
field, gal = count_field(det)
std = cal_dispersion(gal, det)
print('all', std, len(field)/50 * 100, len(det)/len(field))

judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
field, gal = count_field(det)
std = cal_dispersion(gal, det)
print('3031', std, len(field)/11 * 100, len(det)/len(field))

det = pd.read_csv('../detection/ghosts_detection_v4.4.csv')
judge = det['RGB AGB Ratio'] >= 4
det = det[judge].reset_index(drop=True)
field, gal = count_field(det)
std = cal_dispersion(gal, det)
print('all', std, len(field)/50 * 100, len(det)/len(field))

judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
field, gal = count_field(det)
std = cal_dispersion(gal, det)
print('3031', std, len(field)/11 * 100, len(det)/len(field))

det = pd.read_csv('../detection/ghosts_detection_v4.4.csv')
judge = (det['RGB AGB Ratio'] >= 4) & (det['# star below tip'] >= 200)
det = det[judge].reset_index(drop=True)
field, gal = count_field(det)
std = cal_dispersion(gal, det)
print('all', std, len(field)/50 * 100, len(det)/len(field))

judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
field, gal = count_field(det)
std = cal_dispersion(gal, det)
print('3031', std, len(field)/11 * 100, len(det)/len(field))

det = pd.read_csv('../detection/ghosts_detection_v4.1.csv')
judge = (det['RGB AGB Ratio'] >= 4) & (det['# star below tip'] >= 200)
det = det[judge].reset_index(drop=True)
field, gal = count_field(det)
std = cal_dispersion(gal, det)
print('all', std, len(field)/50 * 100, len(det)/len(field))

judge = ['3031' in det['field'][i] for i in range(len(det))]
det = det[judge].reset_index(drop=True)
field, gal = count_field(det)
std = cal_dispersion(gal, det)
print('3031', std, len(field)/11 * 100, len(det)/len(field))
