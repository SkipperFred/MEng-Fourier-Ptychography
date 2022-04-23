# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 15:51:21 2022

@author: freds
"""
import os
import cv2
import numpy as np
import math

def roundHalfUp(num):
    decVal = num-math.floor(num)
    if (decVal<0.5):
        roundNum = math.floor(num)
    elif (decVal>=0.5):
        roundNum = math.ceil(num)
    return roundNum

def LED_Location(xstart,ystart,arraysize):
    xorynode=np.zeros(70)
    xorynode[0]= 1
    dif = 1
    dif_judge=1
    for i in range(1,70):
        xorynode[i]=xorynode[i-1]+dif
        if (dif_judge<2):
            dif_judge=dif_judge+1
        else:
            dif=dif+1
            dif_judge=1
    
    xlocation = np.zeros(arraysize*arraysize)
    ylocation = np.zeros(arraysize*arraysize)
    xlocation[0] = xstart
    ylocation[0] = ystart
    xy_order = 2
    
    for i in range(1,arraysize*arraysize):
        if (i+1<=xorynode[xy_order-1]):
            pass
        else:
            xy_order = xy_order+1
        
        if ((xy_order%2)==0):
            xlocation[i] = xlocation[i-1]+(int(-1)**int(((xy_order/2)%2)+1))
            ylocation[i] = ylocation[i-1]
        elif ((xy_order%2)==1):
            xlocation[i] = xlocation[i-1]
            ylocation[i] = ylocation[i-1]+(int(-1)**int((((xy_order-1)/2)%2)+1))
    return xlocation, ylocation  

instances = []

# Load in the images
for filepath in os.listdir('Experiment DataSet\paper USAF 2104/'):
    instances.append(cv2.imread('Experiment DataSet\paper USAF 2104/'+filepath,1))

dataSet = np.asarray(instances)
greenDataSet = dataSet[:,:,507:3547,1]
greenDataSet = np.transpose(greenDataSet, [1,2,0])

theta = 0
wavelength = 5.32e-07

xint = 0
yint = 0
edge = 60

arraysize = 15
LEDheight      = 60;
LEDgap   = 3.4*0.5
pratio = 2

xlocation = np.zeros([arraysize,arraysize])
ylocation = np.zeros([arraysize,arraysize])

for i in range(arraysize):
    xlocation[i,:] = np.arange(-(arraysize-1)/2,(arraysize)/2, 1)*(LEDgap)
    ylocation[:,i] = np.arange(-(arraysize-1)/2,(arraysize)/2, 1)*(LEDgap)
    
xlocation = np.reshape(xlocation, [1,arraysize**2])
ylocation = np.reshape(ylocation, [1,arraysize**2])

kx_relative = np.sin(np.arctan(xlocation/LEDheight));
ky_relative = np.sin(np.arctan(ylocation/LEDheight));

k0 = 2*np.pi/wavelength;
kx = k0*kx_relative;
ky = k0*ky_relative;
spsize = 1.27e-6;
psize = spsize/pratio;
NA = 0.4;
used_img = np.arange(0,64,1);

m1 =  np.size(greenDataSet,0);
n1 =  np.size(greenDataSet,1);

m = int((spsize/psize)*m1);
n = int((spsize/psize)*n1);

dst = np.zeros(shape=(m,n))

dkx = 2*np.pi/(psize*n);
dky = 2*np.pi/(psize*m);
cutoffFrequency = NA*k0;
kmax = np.pi/spsize;

[kxm,kym]=np.meshgrid(np.arange(-kmax,kmax+1,kmax/((n1-1)/2)),np.arange(-kmax,kmax+1,kmax/((n1-1)/2)));

CTF = np.zeros([m1,n1]);
CTF = ((np.power(kxm,2)+np.power(kym,2))<cutoffFrequency**2);

xstart = 8 
ystart = 8
[xlocation2, ylocation2] = LED_Location(xstart, ystart, arraysize)
xlocation2 = np.subtract(xlocation2,1)
xlocation2 = xlocation2.astype(np.int32)
ylocation2 = np.subtract(ylocation2,1)
ylocation2 = ylocation2.astype(np.int32)
ylocation2 = np.subtract(arraysize-1, ylocation2)
seq = np.zeros([1,arraysize**2], np.int8)

grid = np.arange(0,arraysize**2)
grid = np.reshape(grid,[arraysize, arraysize])


for i in range(0,max(xlocation2.shape)):
    seq[0,i] = int(grid[ylocation2[i], xlocation2[i]])

objectRecover = np.ones([m,n])
objectRecoverFT = np.fft.fftshift(np.fft.fft2(objectRecover));
loop = 5;
for tt in range(loop):
    print(tt)
    for i3 in range(0,max(used_img.shape)):
        print(i3)
        i2 = seq[0,i3]
        kxc = roundHalfUp((n+1)/2 + kx[0,i2]/dkx)
        kyc = roundHalfUp((m+1)/2 + ky[0,i2]/dky)
        kyl = roundHalfUp(kyc-(m1-1)/2)
        kyh = roundHalfUp(kyc+(m1-1)/2)
        kxl = roundHalfUp(kxc-(n1-1)/2)
        kxh = roundHalfUp(kxc+(n1-1)/2)
        
        lowResFT = ((m1/m)**2)*np.multiply(objectRecoverFT[kyl:kyh+1,kxl:kxh+1],CTF);
        # lowResFT = np.zeros([m1,n1], np.complex128)
        # for a in range(kyl, kyh+1):
        #     for b in range(kxl, kxh+1):
        #         lowResFT[(a-kyl),(b-kxl)] = ((m1/m)**2)*np.multiply(objectRecoverFT[a,b],CTF[a-kyl,b-kxl])
        
        im_lowRes = np.fft.ifft2(np.fft.ifftshift(lowResFT));
        im_lowRes = ((m1/m)**2)*np.multiply(greenDataSet[:,:,i3],np.exp(np.multiply(1j,np.angle(im_lowRes))));
        lowResFT=np.multiply(np.fft.fftshift(np.fft.fft2(im_lowRes)),CTF);
        objectRecoverFT[kyl:kyh+1,kxl:kxh+1]=np.multiply((1-CTF),objectRecoverFT[kyl:kyh+1,kxl:kxh+1]) + lowResFT;
        # for a in range(kyl, kyh+1):
        #     for b in range(kxl, kxh+1):
        #         objectRecoverFT[a,b] = np.multiply((1-CTF[a-kyl,b-kxl]),objectRecoverFT[a,b]) +lowResFT[a-kyl,b-kxl]

objectRecover = np.fft.ifft2(np.fft.ifftshift(objectRecoverFT))
himOutAmp = cv2.normalize(np.abs(objectRecover[:,:]), None, 255,0, cv2.NORM_MINMAX, cv2.CV_8UC1)
himOutPha =cv2.normalize(np.angle(objectRecover[:,:]),dst,0,1,cv2.NORM_MINMAX)
himOutReal = cv2.normalize(np.real(objectRecover[:,:]), None, 255,0, cv2.NORM_MINMAX, cv2.CV_8UC1)
cv2.imshow("Amplitude", cv2.resize(himOutAmp,(380*2,380*2)))
cv2.imshow("Phase", cv2.resize(himOutPha,(380*2,380*2)))
cv2.imshow("Combined", cv2.resize(himOutReal,(380*2,380*2)))