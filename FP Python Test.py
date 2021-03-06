# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 18:09:45 2021

@author: freds
"""
import numpy as np
from scipy.io import loadmat

dataSet = loadmat('MouseKidney_green')

wlength = dataSet['wlength']
psize = 1.845e-6/4
n = 4*1520
m = 4*1520 
kmax = np.pi/psize
k0 = 2*np.pi/wlength


kx2 = np.arange(-kmax, kmax+ (kmax/((n-1)/2)), (kmax/((n-1)/2)))
    
kyVal = -kmax
ky2 = []
ky2.append(kyVal)
while kyVal <= kmax:
    kyVal = kyVal+ (kmax/((m-1)/2))
    ky2.append(kyVal)
    
[kxm, kym] = np.meshgrid(kx2,ky2)
kzm = np.sqrt(k0**2-np.power(kxm, 2)-np.power(kym,2))