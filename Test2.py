# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 11:19:09 2022

@author: freds
"""
import numpy as np

def cart2pol(x,y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y,x)
    return (phi,rho)

def pol2cart(phi,rho):
    x = rho*np.cos(phi)
    y = rho*np.sin(phi)
    return(x, y)

m =2
n=2
m1 = 201
n1 = 201
tpixel = max(m1,n1)
NAfily = 69.7077
NAfilx = 69.7077
NApixel = 2*max(round(NAfily),round(NAfilx))

xVal = -tpixel/NApixel
x = []
x.append(xVal)
while xVal <= tpixel/NApixel:
    xVal = xVal+ (((tpixel/NApixel)-(-tpixel/NApixel))/(tpixel-1))
    x.append(xVal)
x = np.reshape(x, (1,tpixel))

[X,Y] = np.meshgrid(x,x)
[theta, r] = cart2pol(X,Y)
idx = np.zeros((tpixel,tpixel),dtype=int)
for a in range(tpixel):
    for b in range(tpixel):
        if (r[a,b]<=1):
            idx[a,b] = 1
        else:
            idx[a,b] = 0
            
z = np.zeros(np.shape(X))
r = np.multiply(r,idx)
r = r[np.nonzero(r)]
r = np.reshape(r,(len(r),1))

theta = np.multiply(theta, idx)
theta = theta[np.nonzero(theta)]
theta = np.reshape(theta,(len(theta),1))

n = [n]
m_abs = [abs(m)]
rpowers = []
length_r = len(r)

for j in range(len(n)):
    if (m_abs[j]-n[j] ==0):
        rpowers.append(2)
    else:
        for i in range(m_abs[j],n[j],2):
            rpowers.append(i)

rpowers = np.unique(rpowers)
rpowern = np.zeros(np.shape(r))

if (rpowers[0]==0):
    for i in range(1, len(rpowers)):
        for a in range(length_r):
            rpowern[a,0] = np.power(r[a,0],rpowers[i])
else:
    for i in range(0, len(rpowers)):
        for a in range(length_r):
            rpowern[a,0] = np.power(r[a,0],rpowers[i])

z = np.zeros((length_r, len(n)))
for j in range(len(n)):
    s = []
    for a in range(int((n[j]-m_abs[j]/2))):
        s.append(a)
    pows=[]    
    if (n[j]-m_abs[j] ==0):
        pows.append(2)
    else:
        for i in range(n[j],m_abs[j],-2):
            pows.append(i)
    if (len(s)-1 ==0):
        p=(1-2*(s[0]%2))*(np.prod(range(2,n[j]-s[0]))+1)/(np.prod(range(2,s[0])))/np.prod(range(2,(int((n[j]-m_abs[j])/2)-s[0])))/(np.prod(range(2,(int((n[j]+m_abs[j])/2)-s[0])))+1)
        idx = (pows[0]==rpowers)
        z[:,j] = z[:,j]+p*rpowern[:,0]*idx
    else:
        for i in range(len(s)-1,0,-1):
           p=(1-2*(s[0]%2))*(np.prod(range(2,n[j]-s[0]))+1)/(np.prod(range(2,s[0])))/np.prod(range(2,(int((n[j]-m_abs[j])/2)-s[0])))/(np.prod(range(2,(int((n[j]+m_abs[j])/2)-s[0])))+1)
           idx = (pows[0]==rpowers)
           z[:,j] = z[:,j]+p*rpowern[:,0]*idx

idx_pos = m>0
idx_neg = m<0
if np.any(idx_pos):
    z[:,idx_pos-1] = abs(np.multiply(z[:,idx_pos-1], (np.asmatrix(np.cos(theta*m_abs[idx_pos-1])).H)))
if np.any(idx_neg):
    z[:,idx_neg-1] = abs(np.multiply(z[:,idx_neg-1], (np.asmatrix(np.cos(theta*m_abs[idx_neg-1])).H)))
                               

