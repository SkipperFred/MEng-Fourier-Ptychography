# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 22:57:09 2021

@author: freds
"""
from scipy.io import loadmat
import cv2
import numpy as np
import math
from imresize import imresize


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

def calculate(x0,y0,H,h,n):
    l = math.sqrt(x0*x0+y0*y0)
    thetal = math.atan2(y0,x0)
    
    xoff = 0;
    thetag = -math.asin(l/math.sqrt(l*l+H*H)/n)
    xint = h*math.tan(thetag)
    xoff = xoff-xint
    
    while abs(xint) > 0.001:
        thetag = -math.asin((l-xoff)/math.sqrt((l-xoff)*(l-xoff)+H*H)/n)
        xint = xoff+h*math.tan(thetag)
        xoff = xoff-xint
    
    theta = math.asin((l-xoff)/math.sqrt((l-xoff)*(l-xoff)+H*H))
    
    NAt = abs(math.sin(theta))
    kx = -NAt*math.cos(thetal)
    ky = -NAt*math.sin(thetal)
    
    return kx, ky, NAt

def k_vector(xi, yi, H, LEDp, nglass, t, theta, xint, yint, total):
    kx = np.zeros((1,total))
    ky = np.zeros((1,total))
    NAt = np.zeros((1,total))
    
    for tt in range(total):
        x0 = xint+xi[0,tt]*LEDp #from rightmost position
        y0 = yint+yi[0,tt]*LEDp #from topmost position
        x1 = x0*math.cos(theta*math.pi/180)-y0*math.sin(theta*math.pi/180)
        y1 = x0*math.sin(theta*math.pi/180)+y0*math.cos(theta*math.pi/180)
        kx[0,tt],ky[0,tt],NAt[0,tt] = calculate(x1,y1,H,t,nglass)
    return kx, ky, NAt

def cart2pol(x,y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y,x)
    return (phi,rho)

def pol2cart(phi,rho):
    x = rho*np.cos(phi)
    y = rho*np.sin(phi)
    return(x, y)

def zernfun(n,m,r,theta):
    
    r = r[np.nonzero(r)]
    r = np.reshape(r,(len(r),1))
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
        rpowern= np.insert(rpowern,2)
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
    return z

def gzn(tpixel, NApixel,m,n):
    
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
    z = np.ones(np.shape(X))
    z = np.multiply(z,idx)
    
    zOut = zernfun(n,m, np.multiply(r,idx), np.multiply(theta, idx))
    currCount = 0
    for a in range(len(z)):
        for b in range(len(z[0])):
            if (z[a,b]!=0):
                z[a,b] = zOut[currCount,0]
                currCount= currCount+1
    return z
        

def himrecover(imseqlow,kx,ky,NA,wlength,spsize,psize,z, opts):
    
    loopNum = opts["loopNum"]
    alpha = opts["alpha"]
    beta = opts["beta"]
    gamma_obj = opts["gamma_obj"]
    gamma_p = opts["gamma_p"]
    eta_obj = opts["eta_obj"]
    eta_p = opts["eta_p"]
    T = opts["T"]
    aberration = opts["aberration"]
          
    #k-space parameterisation
    [m1, n1, numim] = imseqlow.shape
    pratio = round(spsize/psize)
    m = pratio*m1
    n = pratio*n1
    k0 = 2*math.pi/wlength
    kx = k0*kx
    ky = k0*ky
    NAfilx = NA*(1/wlength)*n*psize
    NAfily = NA*(1/wlength)*m*psize 
    kmax = math.pi/psize
    dkx = 2*math.pi/(psize*n)
    dky = 2*math.pi/(psize*m)
    
    kx2 = np.arange(-kmax, kmax, (kmax/((n-1)/2)))
    if (kx2.size < n):
        np.append(kx2, kmax)
    
    ky2 = np.arange(-kmax, kmax, (kmax/((n-1)/2)))
    if (ky2.size < m):
        np.append(ky2, kmax)
    
    [kxm, kym] = np.meshgrid(kx2,ky2)
    kzm = np.sqrt(k0**2-np.power(kxm, 2)-np.power(kym,2))
    
    H2 = np.exp(1j*z*kzm.real)*np.exp(-abs(z)*abs(kzm.imag))
    astigx = 0
    astigy = 0
    
    [M1, N1] = np.meshgrid(list(range(1,m1+1)),list(range(1,n1+1)))
    
    zn = astigx*gzn(max(m1,n1),2*max(round(NAfily[0,0]), round(NAfilx[0,0])),2,2)+astigy*gzn(max(m1,n1),2*max(round(NAfily[0,0]), round(NAfilx[0,0])),-2,2)
    zn = np.resize(zn, (m1,n1))
    
    zn = astigx*gzn(max(m1,n1),2*max(round(NAfily[0,0]), round(NAfilx[0,0])),2,2)+astigy*gzn(max(m1,n1),2*max(round(NAfily[0,0]), round(NAfilx[0,0])),-2,2)
    zn = np.resize(zn, (m1,n1))
    
    if (aberration != 0):
        fmaskpro = aberration
    else:
        fmaskproPT1 = np.multiply(1, (np.less_equal((np.power((N1-(m1+1)/2)/NAfily,2))+(np.power((M1-(n1+1)/2)/NAfilx,2)),1)).astype('float64'))
    
    fmaskA = (np.arange(math.ceil((m+1)/2 - (m1-1)/2), (math.ceil((m+1)/2+(m1-1)/2))+1))
    fmaskB = (np.arange((math.ceil((n+1)/2 - (n1-1)/2)), (math.ceil((n+1)/2+(n1-1)/2))+1))
    fmaskproPT2 = np.zeros((len(fmaskA),len(fmaskB)), dtype="complex128")
    for a in range(len(fmaskA)):
        for b in range(len(fmaskB)):
            fmaskproPT2[a,b] = H2[fmaskA[a]-1,fmaskB[b]-1]
    fmaskproPT3 = np.exp(np.multiply(np.pi*1j,zn))
    fmaskpro = np.multiply(np.multiply(fmaskproPT1,fmaskproPT2), fmaskproPT3)

    him = imresize(np.sum(imseqlow, axis = 2), output_shape=(m,n))
    himFT = np.fft.fftshift(np.fft.fft2(him))
    O_j = np.zeros((m1,n1),dtype=complex)
    
    #main part to optimise estimate of high-res image
    for i in range(1,3):
        for i3 in range(0,numim):
            kxc = int(round((n+1)/2-kx[0,i3]/dkx))
            kyc = int(round((m+1)/2-ky[0,i3]/dky))
            
            kyl = int(round(kyc-(m1-1)/2))
            kyh = int(round(kyc+(m1-1)/2))
            
            kxl = int(round(kxc-(n1-1)/2))
            kxh = int(round(kxc+(n1-1)/2))
            
            for a in range(kyl, kyh+1):
                for b in range(kxl, kxh+1):
                    O_j[(a-kyl),(b-kxl)] = himFT[a,b]
            
            lowFT = np.multiply(O_j,fmaskpro)
            im_lowFT = np.fft.ifft2(np.fft.ifftshift(lowFT))
            updatetemp = (pratio**2)*imseqlow[:,:,i3]
            im_lowFT = np.multiply(updatetemp, np.exp(np.multiply(1j,np.angle(im_lowFT))))
            lowFT_p = np.fft.fftshift(np.fft.fft2(im_lowFT))
            
            ftUp = np.conj(fmaskpro)
            ftUp2 = np.max(np.max(np.power((abs(fmaskpro)), 2)))
            ftUp3 = (lowFT_p - lowFT)
            ftTemp = np.multiply((np.divide(ftUp, ftUp2)), ftUp3)
            for a in range(kyl, kyh+1):
                for b in range(kxl, kxh+1):
                    himFT[a,b] = himFT[a,b] + ftTemp[a-kyl,b-kxl]
             
    countimg = 0
    tt = np.reshape(np.ones(loopNum*numim), (1, loopNum*numim))
    
    vobj0 = np.zeros((m,n))
    vp0 = np.zeros((m1,n1))
    ObjT = himFT
    PT = fmaskpro
    
    for i in range(0, loopNum):
        for i3 in range(0,numim):
            countimg=countimg+1
            kxc=int(round((n+1)/2-kx[0,i3]/dkx));
            kyc=int(round((m+1)/2-ky[0,i3]/dky));
            kyl=int(round(kyc-(m1-1)/2));
            kyh=int(round(kyc+(m1-1)/2));  
            kxl=int(round(kxc-(n1-1)/2));
            kxh=int(round(kxc+(n1-1)/2)); 
            
            for a in range(kyl, kyh+1):
                for b in range(kxl, kxh+1):
                    O_j[(a-kyl),(b-kxl)] = himFT[a,b]
        
            lowFT = np.multiply(O_j, fmaskpro)
            im_lowFT = np.fft.ifft2(np.fft.ifftshift(lowFT))
            tt[0,i3+(i)*numim]=(np.mean(np.mean(abs(im_lowFT)))/np.mean(np.mean((pratio**2)*abs(imseqlow[:,:,i3]))))
        
            if (i>1):
                imseqlow[:,:,i3]=np.multiply(imseqlow[:,:,i3], tt[0,i3+(i)*numim])
        
            updatetemp = np.multiply(pratio**2, imseqlow[:,:,i3])
            im_lowFT = np.multiply(updatetemp, np.exp(np.multiply(1j,np.angle(im_lowFT))))
            lowFT_p = np.fft.fftshift(np.fft.fft2(im_lowFT))
            
            ftUp = (np.multiply((1-alpha),(np.power(abs(fmaskpro),2))))
            ftUp2 = np.multiply(alpha, np.max(np.max(np.power(abs(fmaskpro),2))))
            ftUp3 = (lowFT_p-lowFT)
            ftUp4 = np.conj(fmaskpro)
            ftTemp = np.multiply(gamma_obj, np.multiply(ftUp4,(np.divide(ftUp3, ftUp + ftUp2))))
        
            for a in range(kyl,kyh+1):
                for b in range(kxl, kxh+1):
                    himFT[a,b] = himFT[a,b] + ftTemp[a-kyl,b-kxl]
        
            fmaskUp = np.multiply(beta,np.max(np.max(np.power(abs(O_j),2))))
            fmaskUp2 = np.multiply((1-beta),np.power(abs(O_j),2))
            fmaskUp3 = (lowFT_p-lowFT)
            fmaskUp4 = np.conj(O_j)
            fmaskTemp = np.multiply(gamma_p, np.multiply(fmaskUp4,np.divide(fmaskUp3,(fmaskUp2 + fmaskUp)) ))
        
            for a in range(0, m1):
                for b in range(0,n1):
                    fmaskpro[a,b] = fmaskpro[a,b]+ fmaskTemp[a,b]
            
            if (countimg == T):
                vobj = np.multiply(eta_obj,vobj0) + (himFT - ObjT)
                himFT = ObjT+vobj0
                vobj0 = vobj
                ObjT = himFT   
                vp = np.multiply(eta_p,vp0) + (fmaskpro - PT)
                fmaskpro = PT + vp
                vp0 = vp
                PT = fmaskpro
                
                countimg = 0
    
    him = np.fft.ifft2(np.fft.ifftshift(himFT))
    return him, tt, fmaskpro, imseqlow


dataSet = loadmat('MouseKidney_green')
imgs = dataSet['imlow_HDR']
dst = np.zeros(shape=(201,201))

is_show = 'centre' #'centre' shows first low res raw image; 'all' loads all dynamically
if (is_show == 'centre'):
    cv2.imshow("Mouse Kidney", imgs[:,:,1])
elif(is_show=='all'):
    for i in range(225):
        cv2.imshow("Mouse Kidney",imgs[:,:,i])
        cv2.waitKey()
        cv2.destroyAllWindows()

#Setup experiment Parameters   
xstart = 18 #absolute coordinate of initial LED
ystart = 20
arraysize = 15 # side length of lit LED array
[xlocation, ylocation] = LED_Location(xstart, ystart, arraysize)

xlocation = np.reshape(xlocation, (1,225))
ylocation = np.reshape(ylocation, (1,225))

H = 90.88 #distance between LED and sample in mm
LEDp = 4 #distance between adjacent LEDs, in mm
nglass = 1.52 #refraction index of glass substrate
t = 1 #glass thickness in mm
theta = 0 #rotation angle of LED array to the camera sensor frame, in degrees
xint = 0 # off set of initial LED to the patch centre, in mm
yint = 0

[kx, ky, NAt] = k_vector(xlocation-xstart, ylocation-ystart, H, LEDp, nglass, t, theta, xint, yint, arraysize*arraysize)

#Reconstruct by FP algorithm
NA = 0.1
spsize = 1.845e-6
upsmp_ratio = 4
psize = spsize/upsmp_ratio

wlength = dataSet['wlength']
z = dataSet['z']

opts = {
    "loopNum" : 10, #iteration number
    "alpha" : 1, #'1' for ePIE, other value for rPIE
    "beta" : 1, #'1' for ePIE, other value for rPIE
    "gamma_obj" : 1, #step size for object updating
    "gamma_p" : 1, #step size for pupil updating
    "eta_obj" : 0.2, #stepsize for adding momentumn to object updating
    "eta_p" : 0.2, #stepsize for adding momentumn to pupil updating
    "T" : 1, #do momentumn every T images
    "aberration" : dataSet['aberration'], #pre-calibrated aberration, if available
}





used_idx = list(range(0,arraysize**2))


imlow_used = imgs[:,:,used_idx]

kx_used = kx[0,used_idx]
ky_used = ky[0,used_idx]

[him, tt, fprobe, imLow_HDR1] = himrecover(imlow_used, kx_used, ky_used, NA, wlength, spsize, psize, z, opts)

himOutAmp =cv2.normalize(np.abs(him[:,:]),dst,0,1,cv2.NORM_MINMAX)
himOutPha =cv2.normalize(np.angle(him[:,:]),dst,0,1,cv2.NORM_MINMAX)
cv2.imshow("Mouse Kidney Amplitude", himOutAmp)
cv2.imshow("Mouse Kidney Phase", himOutPha)
    