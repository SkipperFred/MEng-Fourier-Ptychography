# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 21:42:51 2022

@author: freds
"""
import numpy as np
from sender import Sender
import time

arraysize = 8
xstart = 3
ystart = 3

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

xlocation = np.zeros(arraysize*arraysize, dtype=int)
ylocation = np.zeros(arraysize*arraysize, dtype=int)
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
        
lookUpTable = list(range(0,arraysize**2))
lookUpTable = np.reshape(lookUpTable, (arraysize, arraysize))

s = Sender('COM3')


for i in range(0,arraysize**2):
    print(i)
    LEDnum = lookUpTable[xlocation[i],ylocation[i]]
    s.send('LEDOn('+str(LEDnum)+')')
    time.sleep(1)
    s.send('LEDOff('+str(LEDnum)+')')
    
    
s.close()