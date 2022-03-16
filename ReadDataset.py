# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 14:14:16 2022

@author: freds
"""

import numpy as np
import cv2
import os

def cropImage (img, cropfactor):
    num, cropy, cropx = img.shape
    cropystart = int(cropy/cropfactor - cropy/(cropfactor*2))
    cropyend = int(cropy/cropfactor + cropy/(cropfactor*2))
    cropxstart = int(cropx/cropfactor - cropy/(cropfactor*2))
    cropxend = int(cropx/cropfactor + cropy/(cropfactor*2))
    
    img = img[:, cropystart:cropyend, cropxstart:cropxend]
    return img

instances = []

# Load in the images
for filepath in os.listdir('Experiment DataSet\img_database_jpg_2022_03_09_1/'):
    instances.append(cv2.imread('Experiment DataSet\img_database_jpg_2022_03_09_1/'+filepath,1))

dataSet = np.asarray(instances)
greenDataSet = dataSet[:,:,:,1]

greenDataSet = cropImage(greenDataSet, 2)


cv2.imshow("Data", cv2.resize(greenDataSet[0,:,:],(380,380)))