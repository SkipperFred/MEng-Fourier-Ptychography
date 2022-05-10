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

def calculateRing(arrayAngleOffSet,LEDNum,H):
    arrayAngleOffSet = arrayAngleOffSet*np.pi/180
    if LEDNum<12:
        l = 21
        thetal = (12-LEDNum)*(2*np.pi/12)+arrayAngleOffSet
    elif LEDNum<12+16:
        l = 30
        thetal = (16-(LEDNum-12))*(2*np.pi/16)+arrayAngleOffSet
    elif LEDNum<12+16+24:
        l = 39
        thetal = (24-(LEDNum-(12+16)))*(2*np.pi/24)+arrayAngleOffSet
    elif LEDNum<12+16+24+32:
        l = 48
        thetal = (32-(LEDNum-(12+16+24)))*(2*np.pi/32)+arrayAngleOffSet
    
    theta = math.atan2(l, H)

    NAt = abs(math.sin(theta))
    kx = l*np.sin(thetal)
    ky = l*np.cos(thetal)
    
    kx_relative = np.sin(np.arctan(kx/LEDheight));
    ky_relative = np.sin(np.arctan(ky/LEDheight));
    
    k0 = 2*np.pi/wavelength;
    kx = k0*kx_relative;
    ky = k0*ky_relative;
    
    return kx, ky, NAt

def k_vectorRing(xi, H, xint, yint, total, angleOffset):
    kx = np.zeros((1,total))
    ky = np.zeros((1,total))
    NAt = np.zeros((1,total))

    for tt in range(total):
        kx[0,tt],ky[0,tt],NAt[0,tt] = calculateRing(angleOffset,tt,H)    
    return kx, ky, NAt

def loadImage(filepathName):
    instances = []

    # Load in the images
    for filepath in os.listdir(filepathName):
        instances.append(cv2.imread(filepathName+filepath,1))#Images loaded as numpy arrays into list

    dataSet = np.asarray(instances) #Convert list into numpy array

    #find the colour channel with highest mean value
    if (np.mean(dataSet[:,:,:,0])>np.mean(dataSet[:,:,:,1])) and (np.mean(dataSet[:,:,:,0])>np.mean(dataSet[:,:,:,2])):
        colourChannel = 0
        print("r")
    elif (np.mean(dataSet[:,:,:,1])>np.mean(dataSet[:,:,:,0])) and (np.mean(dataSet[:,:,:,1])>np.mean(dataSet[:,:,:,2])):
        colourChannel = 1;
        print("g")
    elif (np.mean(dataSet[:,:,:,2])>np.mean(dataSet[:,:,:,0])) and (np.mean(dataSet[:,:,:,2])>np.mean(dataSet[:,:,:,1])):
        colourChannel = 2;
        print("b")

    DataSetOut = dataSet[:,:,507:3547,colourChannel] #Sekect only one colour channel
    DataSetOut = np.transpose(DataSetOut, [1,2,0])
    return DataSetOut

DataSet = loadImage('Experiment DataSet\MicroscopeRing2704/')
theta = 0
wavelength = 5.32e-07

edge = 60

LEDheight = 70
xint = 0
yint = 0
arraysize = 84
wavelength = 5.32e-07
spsize =0.6e-6
psize = spsize/2
angleOffset = -120
NA = 0.5;
used_img = np.arange(0,84-32-24,1);
k0 = 2*np.pi/wavelength;

m1 =  np.size(DataSet,0);
n1 =  np.size(DataSet,1);

m = int((spsize/psize)*m1);
n = int((spsize/psize)*n1);

dst = np.zeros(shape=(m,n))

dkx = 2*np.pi/(psize*n);
dky = 2*np.pi/(psize*m);
cutoffFrequency = NA*k0;
kmax = np.pi/spsize;

[kxm,kym]=np.meshgrid(np.arange(-kmax,kmax+1,kmax/((n1-1)/2)),np.arange(-kmax,kmax+1,kmax/((n1-1)/2)));

CTF = np.zeros([m1,n1])
CTF = ((np.power(kxm,2)+np.power(kym,2))<cutoffFrequency**2)


xstart = 8 
ystart = 8
seq = np.asarray(range(arraysize))

[kx, ky, NAt] = k_vectorRing(seq, LEDheight, xint, yint, arraysize, angleOffset)




grid = np.arange(0,arraysize**2)
grid = np.reshape(grid,[arraysize, arraysize])


objectRecover = np.ones([m,n])
objectRecoverFT = np.fft.fftshift(np.fft.fft2(objectRecover));
loop = 5;
for tt in range(loop):
    print(tt)
    for i3 in range(0,max(used_img.shape)):
        print(i3)
        i2 = seq[i3]
        kxc = roundHalfUp((n+1)/2 + kx[0,i2]/dkx)
        kyc = roundHalfUp((m+1)/2 + ky[0,i2]/dky)
        kyl = roundHalfUp(kyc-(m1-1)/2)
        kyh = roundHalfUp(kyc+(m1-1)/2)
        kxl = roundHalfUp(kxc-(n1-1)/2)
        kxh = roundHalfUp(kxc+(n1-1)/2)
        
        lowResFT = ((m1/m)**2)*np.multiply(objectRecoverFT[kyl:kyh+1,kxl:kxh+1],CTF);
        im_lowRes = np.fft.ifft2(np.fft.ifftshift(lowResFT));
        im_lowRes = ((m1/m)**2)*np.multiply(DataSet[:,:,i3],np.exp(np.multiply(1j,np.angle(im_lowRes))));
        lowResFT=np.multiply(np.fft.fftshift(np.fft.fft2(im_lowRes)),CTF);
        objectRecoverFT[kyl:kyh+1,kxl:kxh+1]=np.multiply((1-CTF),objectRecoverFT[kyl:kyh+1,kxl:kxh+1])+lowResFT
       
objectRecover = np.fft.ifft2(np.fft.ifftshift(objectRecoverFT))
himOutAmp = cv2.normalize(np.abs(objectRecover[:,:]), dst, 255,0, cv2.NORM_MINMAX, cv2.CV_8UC1)
himOutPha =cv2.normalize(np.angle(objectRecover[:,:]),dst,0,1,cv2.NORM_MINMAX)
himOutReal = cv2.normalize(np.real(objectRecover[:,:]), dst, 255,0, cv2.NORM_MINMAX, cv2.CV_8UC1)
logOut = cv2.normalize(np.float32(np.log(objectRecoverFT)),dst, 0,1,cv2.NORM_MINMAX)
cv2.imshow("Amplitude", cv2.resize(himOutAmp,(380*2,380*2)))
cv2.imshow("Phase", cv2.resize(himOutPha,(380*2,380*2)))
cv2.imshow("Combined", cv2.resize(himOutReal,(380*2,380*2)))
cv2.imshow("Log Spectrum", cv2.resize(logOut,(380*2,380*2)))