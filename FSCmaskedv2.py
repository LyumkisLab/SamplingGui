#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 00:46:20 2024

@author: pbaldwin
"""



#%%


# os.chdir('/home/pbaldwin/Dropbox/NYSBC/FSC/Venkat')
#from  EMAN2 import *
import time
import shutil;
import numpy as np;
import matplotlib.pyplot as plt
import numpy
import sys
from sys import argv
import csv
import sys, os
import pandas as pd
import math
from math import *
import matplotlib
import scipy
from scipy import special
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import glob;
import mrcfile


#%%  Section -1 Function Definitions

#@autojit
def ZeroPad(nx,ny,nz,fT,gT):
    fp=np.zeros([nx+2,ny,nz]);
    gp=np.zeros([nx+2,ny,nz]);

    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                fp[ix,iy,iz]=fT[ix,iy,iz]
                gp[ix,iy,iz]=gT[ix,iy,iz]
    return [fp,gp]

#Functions read in 0.078624 seconds for size nx=256 
#%%  Section -1 Function Definitions

#@autojit
def FFTArray2Real(nx,ny,nz,F):

    d1=np.zeros([nx+4,ny,nz]);# 36,32,32

    for ix in range(0,nx+3,2):
        ixover2= ix//2;
        for iy in range(ny):
            for iz in range(nz):
                FNow=F[ixover2,iy,iz];
                d1[ix  ][iy][iz]=FNow.real;
                d1[ix+1][iy][iz]=FNow.imag;
    return d1

#%%  Section -1 Function Definitions

#@autojit
def FFTArray2RealB(nx,ny,nz,F):

    d1=np.zeros([nx+4,ny,nz]);# 36,32,32

    for ix in range(0,nx+3,2):
        ixover2= ix//2;
        for iy in range(ny):
                FNowVec=F[ixover2,iy,:];
                d1[ix  ][iy][:]=FNowVec.real;
                d1[ix+1][iy][:]=FNowVec.imag;
    return d1




#%%


def CreateFSCOutputs(inc,nx,ny,nz,d1,d2):

    # We want the unit vectors corresponding to angles Y
    # We are going to rotate zhat by rotation around y
    # VecUnits
    
    ret = np.zeros(inc+1)
    n1  = np.zeros(inc+1)
    n2  = np.zeros(inc+1)
    lr = np.zeros(inc+1)
    #Count=0;
    for iz in range(nz):
        kz=iz;
        if (iz>nz2): kz=iz-nz;# This is the actual value kz, can be neg
        argz= float(kz*kz)*dz2;
        for iy in range(ny):
            ky=iy;        # This is the actual value ky, can be neg
            if (iy>ny2): ky=iy-ny;
            argy= argz+float(ky*ky)*dy2;
            for ix in range(0,nx,2):
                kx = ix//2;
                # 
                if ( (ix==0) & (kz<0)& (ky<0)):  continue;
                argx = 0.5*np.sqrt(argy + float(ix*ix)*0.25*dx2);
                r=int(round(inc*2*argx));
                if (r <= inc):
                    #ii = ix + (iy  + iz * ny)* lsd2;
                    ret[r] += d1[ix  ,iy,iz] * d2[ix,iy,iz] ;
                    ret[r] += d1[ix+1,iy,iz] * d2[ix+1,iy,iz] ;
                    n1[r]  += d1[ix  ,iy,iz] * d1[ix,iy,iz];
                    n1[r]  += d1[ix+1,iy,iz] * d1[ix+1,iy,iz];
                    n2[r]  += d2[ix  ,iy,iz] * d2[ix,iy,iz];
                    n2[r]  += d2[ix+1,iy,iz] * d2[ix+1,iy,iz];
                    lr[r]  += 2.0;# Number of Values

    return [ret,n1,n2,lr]


#%%

# anglesY= np.linspace(-60,60,61, endpoint=True)

def CreateAllFSCOutputs(inc,nx,ny,nz,d1,d2,anglesY):

    # We want the unit vectors corresponding to angles Y
    # We are going to rotate zhat by rotation around y
    
    
    NumAngles = len(anglesY)
    #VecUnits = np.zeros([NumAngles,2])
    VecUnitsCrosszhat = np.zeros([NumAngles,2])
    
    for jAngle in range(NumAngles):
        AngleNow = anglesY[jAngle] *np.pi/180.0
        #VecUnits[jAngle,:] =[np.cos(AngleNow),np.sin(AngleNow) ] 
        VecUnitsCrosszhat[jAngle,:] =[np.sin(AngleNow),-np.cos(AngleNow) ] 
        
    # VecUnits
    # we want Vecunits
    
    retIn = np.zeros(inc+1) ;    retOut = np.zeros(inc+1) ;
    n1In  = np.zeros(inc+1);     n1Out = np.zeros(inc+1) ;
    n2In  = np.zeros(inc+1);     n2Out = np.zeros(inc+1) ;
    lrIn  = np.zeros(inc+1);     lrOut = np.zeros(inc+1) ;
    #Count=0;
    for iz in range(nz):
        kz=iz;
        if (iz>nz2): kz=iz-nz;# This is the actual value kz, can be neg
        argz= float(kz*kz)*dz2;
        for iy in range(ny):
            ky=iy;        # This is the actual value ky, can be neg
            if (iy>ny2): ky=iy-ny;
            argy= argz+float(ky*ky)*dy2;
            for ix in range(0,nx,2):
                kx = ix//2;
                # 
                if ( (ix==0) & (kz<0)& (ky<0)):  continue;
                argx = 0.5*np.sqrt(argy + float(ix*ix)*0.25*dx2);
                r=int(round(inc*2*argx));
                NHat = np.array([kx,kz])
                if (r <= inc):
                    Distances = np.dot(VecUnitsCrosszhat,NHat)
                    vv=np.where(np.abs(Distances) <0.5)[0]
                    
                    if len(vv):
                        #ii = ix + (iy  + iz * ny)* lsd2;
                        retIn[r] += d1[ix  ,iy,iz] * d2[ix,iy,iz] ;
                        retIn[r] += d1[ix+1,iy,iz] * d2[ix+1,iy,iz] ;
                        n1In[r]  += d1[ix  ,iy,iz] * d1[ix,iy,iz];
                        n1In[r]  += d1[ix+1,iy,iz] * d1[ix+1,iy,iz];
                        n2In[r]  += d2[ix  ,iy,iz] * d2[ix,iy,iz];
                        n2In[r]  += d2[ix+1,iy,iz] * d2[ix+1,iy,iz];
                        lrIn[r]  += 2.0;# Number of Values
                    else:
                        retOut[r] += d1[ix  ,iy,iz] * d2[ix,iy,iz] ;
                        retOut[r] += d1[ix+1,iy,iz] * d2[ix+1,iy,iz] ;
                        n1Out[r]  += d1[ix  ,iy,iz] * d1[ix,iy,iz];
                        n1Out[r]  += d1[ix+1,iy,iz] * d1[ix+1,iy,iz];
                        n2Out[r]  += d2[ix  ,iy,iz] * d2[ix,iy,iz];
                        n2Out[r]  += d2[ix+1,iy,iz] * d2[ix+1,iy,iz];
                        lrOut[r]  += 2.0;# Number of Values
                        

    return [retIn,n1In,n2In,lrIn,retOut,n1Out,n2Out,lrOut]


#%%
#def CalculateFSCWithMask(AngleFiley, FN1, FN2,KmaxInPixels):

def GetFFTs(FN1, FN2):
    
    h5f_HalfMap1= mrcfile.open(FN1)
    f= h5f_HalfMap1.data.T
    
    h5f_HalfMap2= mrcfile.open(FN2)
    g= h5f_HalfMap2.data.T
    
    h5f_HalfMap1.close()
    h5f_HalfMap2.close()
    
    
        
    startTime = time.time()
    
    #fT=f.T;# Now it is like EMAN
    #gT=g.T;#
    
    
    
    #[nx,ny,nz] =fT.shape
    [nx,ny,nz] =f.shape
    
    print(nx,ny,nz)
    
    deltaTime =time.time()-startTime;
    print("Maps read in %f seconds for size nx=%g,ny=%g, nz=%g " % (deltaTime,nx,ny,nz))
    


    #[fp,gp]=ZeroPad(nx,ny,nz,fT,gT)
    [fp,gp]=ZeroPad(nx,ny,nz,f,g)





    startTime = time.time()

    F =np.fft.fftn(fp);#F= temp.transpose();
    G =np.fft.fftn(gp);#G= temp.transpose();

    H=F*np.conj(G);
    print(H.shape)
    #HAngle= np.angle(H);
    #CosHangle=np.cos(HAngle);


    deltaTime =time.time()-startTime;
    print("FFTs performed in %f seconds for size nx=%g " % (deltaTime,nx))
    
    return F,G,H,nx,ny,nz

    
    
#%%





print('Hello section 0: Define Variables')

#sys.exit()


#
# 
# FSCMasked.py
#                 HalfMap1.mrc   HalfMap2.mrc  OutputStringLabel   APixels AngleList                      
#
# Creates   ResEMOutresultAve+OutputStringLabel.csv (which is the usual FSC)
#           ResEMOut+OutputStringLabel.hdf which is the 3D FSC file ResEMR (the real part of the cccs) 
#           Plots+OutputStringLabel.jpg which is the slices along x, y, z 


fNHalfMap1='run_half1_class001_unfil.mrc'
fNHalfMap2='run_half2_class001_unfil.mrc'
fNHalfMap1='generate3d_volume_46_41_1_map1.mrc'
fNHalfMap2='generate3d_volume_46_41_1_map2.mrc'
fNHalfMap1='volume_46_43_1_map1Sh.mrc'
fNHalfMap2='volume_46_43_1_map2Sh.mrc'

OutputStringLabel='ApoFullMRC'
OutputStringLabel='GSIRD'
OutputStringLabel='SalkOutF'
v = 1.06049
APixels = 2.14
APixels = 1.048
APixels = 1.048


fNHalfMap1=argv[1];
fNHalfMap2=argv[2];

OutputStringLabel= argv[3]; 

APixels = float(argv[4]);

AngleList = np.loadtxt(argv[5]);


ResultsDir= 'Results'+OutputStringLabel+'/';
isResultsDir=  os.path.isdir(ResultsDir)
if (1-isResultsDir):
    os.mkdir(ResultsDir)

ResEMOutHDF_FN= ResultsDir+'ResEM'+OutputStringLabel+'Out.mrc'
#FTOutHDF_FN= ResultsDir+'FTEM'+OutputStringLabel+'Out.mrc'


if 0:
    dthetaInDegrees =20;
    dthetaInRadians = dthetaInDegrees*np.pi/180.0;
    Thresh = np.cos(dthetaInRadians)   ;
    fractionOfTheSphereAveraged = (1-Thresh)/2;
    # Now, deltaTheta takes up a cone which has area
    #    2 pi * (1-cos(deltaTheta))


    FTOut =  ResultsDir+'FTPlot'+OutputStringLabel;
    PlotsOut= ResultsDir+'Plots'+OutputStringLabel;
    fNResRoot=ResEMOutHDF_FN[0:-4]
    fN_csv_out= fNResRoot+'.csv';

    ResNumOutMRC= fNResRoot+'Num.mrc';# The numerator which is the cross product
    ResDenOutMRC= fNResRoot+'Den.mrc';# The fn for denominator, which is normalization
    
    resultAveOut = fNResRoot+'globalFSC.csv';




FN1= fNHalfMap1
FN2= fNHalfMap2


#%%
#os.chdir('/home/pbaldwin/Downloads/MissingWedgeMikeSchmid/tomograms/')

if 0:
    OutputStringLabel = 'Section7'
    AngleList = np.loadtxt('ListOfAngles.txt')
    FN1 = 'p22_sim_1_gt.mrc'
    FN2 = 'p22_sim_1_eman_new.mrc'; # done  Tu 01/30.2024
    ResultsDir= 'Results'+OutputStringLabel+'/';
    APixels=4.5

# #sys.exit()

# FN1 = 'sph_gt.mrc'
# FN2 = 'sph_imod_resampled.mrc'
# FN2 = 'sph_isonet_nomask_not_resampled.mrc'
# FN2 = 'sph_cn_231128_resampled_no_crowther_scaled_tiltaxis.mrc'


# #if 0: 
#     # FSCMasked.py  HalfMap1.hdf    HalfMap2.hdf   OutputStringLabel   APixels AngleList 
#     #export PATH="/home/software/anaconda3/bin:$PATH"
#     #/home/pbaldwin/Dropbox/BCM2020/IsoNet/FSCmasked.py  sph_gt.mrc   sph_imod_resampled.mrc   Test_FSCmasked   1.0  ListOfAngles.txt 
    
    
    
    
# #AngleFiley=1
# FN1 = 'p22_sim_1_gt.mrc'
# FN2 = 'p22_sim_1_eman_new.mrc'; # done  Tu 01/30.2024
# FN2 = 'p22_sim_1_isonet_nomask_iter15.mrc'; # done  Tu 01/30/2024
# FN2 = 'p2s_1_cn_231128_resampled.mrc'; # done  Tu 01/30/2024
# FN2 = 'p2s_1_cn_231128_lp18-c31_resampled.mrc'; # done  Tu 01/30/2024
# FN2 = 'p22_sim_1_imod_rec_resampled.mrc' ;  # done  Tu 01/30/2024
# #
# #-rw-r--r-- 1 pbaldwin pbaldwin  19376557 Nov 30 04:21 p22_sim_1_eman_new.hdf
# #-rw-r--r-- 1 pbaldwin pbaldwin 186638816 Nov 27 03:00 p22_sim_1_gt.hdf
# #-rw-r--r-- 1 pbaldwin pbaldwin 186625024 Nov 28 13:28 p22_sim_1_imod_rec_resampled.mrc
# #-rw-r--r-- 1 pbaldwin pbaldwin 186625024 Jan 24 13:00 p22_sim_1_isonet_nomask_iter15.mrc
# #-rw-r--r-- 1 pbaldwin pbaldwin  31894696 Nov 21 02:27 p22_sim_1_prjs.hdf



#%%
if 0:
    G_Crowther = G.copy()
  
#%%  Section 1 Perform FFTs
# F,G,H,nz,ny,nz = CalculateFSCWithMask(AngleFiley, FN1, FN2,KmaxInPixels)


print('Hello section 1: Do FFTs')

F,G,H,nx,ny,nz =  GetFFTs(FN1,FN2)



#%%

if 0:
    print(F.shape,G.shape,H.shape)
    HAngle= np.angle(H);
    CosHangle=np.cos(HAngle);

#%%       Section 2 Create Real Arrays as in original EMAN program



print('Hello section 2: Create Real Arrays as in original EMAN program')

startTime = time.time()

d1      = FFTArray2RealB(nx,ny,nz,F);
d2      = FFTArray2RealB(nx,ny,nz,G);
#dcH     = FFTArray2Real(nx,ny,nz,CosHangle);
#dFPower = FFTArray2Real(nx,ny,nz,F*np.conj(F));
#dPR = FFTArray2Real(nx,ny,nz,HAngle);

deltaTime =time.time()-startTime;
print("FFTArray2Real performed in %f seconds for size nx=%g " % (deltaTime,nx))

# FFTArray2Real performed in 22.700670 seconds for size nx=256 , wo autojit
# FFTArray2Real performed in 0.9       seconds for size nx=256 , w autojit
# 

# d1[15][16][17]  is d1.get_value_at(15,16,17)
# But d1[15][16][17], for EMData objects is complex...
# f[16,17,18] = f[16][17][18]
# But this is EMAN's 18,17,16


    

#%%     Section 3  Create FSCs


print('Hello section 3: Create FSC outputs')



nx2 = nx//2;
ny2 = ny//2;
nz2 = nz//2;

anglesY= np.linspace(-60,60,61, endpoint=True)
anglesY = AngleList.copy();



dx2 = 1.0/float(nx2)/float(nx2);
dy2 = 1.0/float(ny2)/float(ny2);
dz2 = 1.0/float(nz2)/float(nz2);
# int inc = Util::round(float(std::max(std::max(nx2,ny2),nz2))/w);
w=1;

inc = max(nx2,ny2,nz2)/w;
inc = int(inc)


startTime = time.time()


#[retcHGlobal,lr]     = CreateFTLikeOutputs(inc,nx,ny,nz,dcH)
[ret,n1,n2,lr] = CreateFSCOutputs(inc,nx,ny,nz,d1,d2);
[retIn,n1In,n2In,lrIn,retOut,n1Out,n2Out,lrOut] = CreateAllFSCOutputs(inc,nx,ny,nz,d1,d2,anglesY);
#[FPower,lr]    = CreateFTLikeOutputs(inc,nx,ny,nz,dFPower)


deltaTime =time.time()-startTime;
print("CreateFSCOutputs performed in %f seconds for size nx=%g " % (deltaTime,nx))

# CreateFSCOutputs performed in  0.393  seconds for size nx=256 ; using autojit
# CreateFSCOutputs performed in 45.944  seconds for size nx=256 ; without autojit
# list(ret)


#%%     Section 4 Create FSCs. Define RMax based on this


print('Hello section 4; create FSCs')


# for nx =32, linc is 17

startTime = time.time()



linc = 0;
for i in range(inc+1):
    if (lr[i]>0):
        linc +=1;

result  = [0 for i in range(3*linc)];
resultIn  = [0 for i in range(3*linc)];
resultOut  = [0 for i in range(3*linc)];


ii = -1;
for i in range(inc+1):
    if (lr[i]>0): 
        ii +=1;
        result[ii]        = float(i)/float(2*inc);
        result[ii+linc]   = float(ret[i] / (np.sqrt(n1[i] * n2[i])));
        result[ii+2*linc] = lr[i]  ;# Number of Values

        resultIn[ii]        = float(i)/float(2*inc);
        resultIn[ii+linc]   = float(retIn[i] / (np.sqrt(n1In[i] * n2In[i])));
        resultIn[ii+2*linc] = lrIn[i]  ;# Number of Values

        resultOut[ii]        = float(i)/float(2*inc);
        resultOut[ii+linc]   = float(retOut[i] / (np.sqrt(n1Out[i] * n2Out[i])));
        resultOut[ii+2*linc] = lrOut[i]  ;# Number of Values


NormalizedFreq = result[0:(inc+1)];
resultAve= result[(inc+1):(2*(inc+1))];# This takes values inc+1 values from inc+1 to 2*inc+1
resultInAve= resultIn[(inc+1):(2*(inc+1))];# This takes values inc+1 values from inc+1 to 2*inc+1
resultOutAve= resultOut[(inc+1):(2*(inc+1))];# This takes values inc+1 values from inc+1 to 2*inc+1


#  filter.lowpass.gauss             :  apix(FLOAT)  
# e2proc3d.py threed.hdf threed.filt.hdf --process=filter.lowpass.gauss:cutoff_freq=0.1
# e2proc3d.py p22_sim_1_gt.mrc p22_sim_1_gt.filt.mrc --process=filter.lowpass.gauss:cutoff_freq=0.01
# e2proc3d.py p22_sim_1_gt.mrc OutputFSC.txt --calcfsc=p22_sim_1_gt.filt.mrc


FSCnow = np.loadtxt('OutputFSC.txt')
plt.plot(FSCnow[:,0],FSCnow[:,1])


deltaTime =time.time()-startTime;
print("Section 4: Writing out FSCs performed in %f seconds for size nx=%g " % (deltaTime,nx))



#%%  Section 5  Make Plots and save to file


print('Hello section 5; Make Plots, save to file')


FN1_Clip = FN1[0:-4]
FN2_Clip = FN2[0:-4]


FN_FSC_All_Measured_Region   = 'FSC_All_'  + FN1_Clip+'_'+FN2_Clip+ '.txt'
FN_FSC_In_Measured_Region    = 'FSC_InMeas_'  + FN1_Clip+'_'+FN2_Clip+ '.txt'
FN_FSC_Outside_Measured_Region   = 'FSC_NotMeas_' + FN1_Clip+'_'+FN2_Clip+ '.txt'
FN_FIg_All                   = 'FSC_All_'+FN1_Clip+'_'+FN2_Clip+ '.jpg'


FN_FSC_All_Measured_Region      =  ResultsDir +  FN_FSC_All_Measured_Region
FN_FSC_In_Measured_Region       =  ResultsDir +  FN_FSC_In_Measured_Region
FN_FSC_Outside_Measured_Region  =  ResultsDir +  FN_FSC_Outside_Measured_Region
FN_FIg_All                      =  ResultsDir +  FN_FIg_All


plt.plot(resultAve)
#plt.title('Usual FSC: Isonet and GT')
plt.plot(resultInAve)
#plt.title('FSC inside wedge: Isonet and GT')
plt.plot(resultOutAve)
plt.title('FSC ')

plt.legend(['Usual FSC','FSC in Measured Region','FSC in Not Measured Region'])

plt.xlim(0,100)
plt.ylim(0,1.05)

#print(FN_FIg_All)
if 0:
    plt.savefig(ResultsDir+'/temp.jpg')
#plt.savefig(FN_FIg_All)




#%%             Section 6 ; write out to file

#os.mkdir(ResultsDir)

print('Hello section 6; write data out to file')


resultAve2File = np.zeros([inc+1,3])
result_In_Measured_2File = np.zeros([inc+1,3])
result_Outside_Measured_2File = np.zeros([inc+1,3])

NormalizedFreq = result[0:(inc+1)];
#APixels=1

for j in range(inc+1):
    valFreqNormalized = NormalizedFreq[j];
    valFreq           = valFreqNormalized/APixels;
    valFSCshell= resultAve[j];
    #AveWriter.writerow([valFreqNormalized,valFreq,valFSCshell])
    resultAve2File[j,:] = [valFreqNormalized,valFreq,valFSCshell]

    valFSCshell= resultInAve[j];
    result_In_Measured_2File[j,:] = [valFreqNormalized,valFreq,valFSCshell]

    valFSCshell= resultOutAve[j];
    result_Outside_Measured_2File[j,:] = [valFreqNormalized,valFreq,valFSCshell]

#k0  = np.array(range(Nxf))/APixels/2.0/Nxf;



np.savetxt(FN_FSC_All_Measured_Region,resultAve2File,fmt='%2.5f')
np.savetxt(FN_FSC_In_Measured_Region,result_In_Measured_2File,fmt='%2.5f')
np.savetxt(FN_FSC_Outside_Measured_Region,result_Outside_Measured_2File,fmt='%2.5f')



print('Exiting from Section 6')



#%%   

print('Hello section 7; replot from files written to disk')


FSC_All_Now = np.loadtxt(FN_FSC_All_Measured_Region)
FSC_Meas_Now = np.loadtxt(FN_FSC_In_Measured_Region)
FSC_UnMeas_Now = np.loadtxt(FN_FSC_Outside_Measured_Region)



plt.figure()
plt.plot(FSC_All_Now[:,1],FSC_All_Now[:,2])
plt.plot(FSC_Meas_Now[:,1],FSC_Meas_Now[:,2])
plt.plot(FSC_UnMeas_Now[:,1],FSC_UnMeas_Now[:,2])

plt.title('FSC vs spatial freq')
plt.legend(['Usual FSC','FSC in Measured Region','FSC in Not Measured Region'])
#if 0:
#    plt.xlim(0,0.2); 
plt.ylim(0,1.1)
plt.xlabel('Inverse Angstroms')
plt.ylabel('FSC')
plt.savefig(FN_FIg_All)




print('Exiting from Section 7')

sys.exit()

#%%

# -rw-r--r-- 1 pbaldwin pbaldwin 13575 Jan 30 13:14 FSC_All_p22_sim_1_gt_p22_sim_1_eman_new.txt
# -rw-r--r-- 1 pbaldwin pbaldwin 13554 Jan 30 13:14 FSC_Out_p22_sim_1_gt_p22_sim_1_eman_new.txt
# -rw-r--r-- 1 pbaldwin pbaldwin 13575 Jan 30 13:14 FSC_In_p22_sim_1_gt_p22_sim_1_eman_new.txt
# -rw-r--r-- 1 pbaldwin pbaldwin 13554 Jan 30 13:31 FSC_Out_p22_sim_1_gt_p22_sim_1_isonet_nomask_iter15.txt
# -rw-r--r-- 1 pbaldwin pbaldwin 13576 Jan 30 13:31 FSC_In_p22_sim_1_gt_p22_sim_1_isonet_nomask_iter15.txt
# -rw-r--r-- 1 pbaldwin pbaldwin 13576 Jan 30 13:31 FSC_All_p22_sim_1_gt_p22_sim_1_isonet_nomask_iter15.txt
# -rw-r--r-- 1 pbaldwin pbaldwin 13575 Jan 30 22:09 FSC_All_p22_sim_1_gt_p2s_1_cn_231128_lp18-c31_resampled.txt
# -rw-r--r-- 1 pbaldwin pbaldwin 13554 Jan 30 22:09 FSC_Out_p22_sim_1_gt_p2s_1_cn_231128_lp18-c31_resampled.txt
# -rw-r--r-- 1 pbaldwin pbaldwin 13575 Jan 30 22:09 FSC_In_p22_sim_1_gt_p2s_1_cn_231128_lp18-c31_resampled.txt
# -rw-r--r-- 1 pbaldwin pbaldwin 13575 Jan 30 22:32 FSC_In_p22_sim_1_gt_p22_sim_1_imod_rec_resampled.txt
# -rw-r--r-- 1 pbaldwin pbaldwin 13575 Jan 30 22:32 FSC_All_p22_sim_1_gt_p22_sim_1_imod_rec_resampled.txt
# -rw-r--r-- 1 pbaldwin pbaldwin 13572 Jan 30 22:32 FSC_Out_p22_sim_1_gt_p22_sim_1_imod_rec_resampled.txt


# IMOD

FSC_All_Now = np.loadtxt('FSC_All_p22_sim_1_gt_p22_sim_1_imod_rec_resampled.txt')
FSC_Meas_Now = np.loadtxt('FSC_In_p22_sim_1_gt_p22_sim_1_imod_rec_resampled.txt')
FSC_UnMeas_Now = np.loadtxt('FSC_Out_p22_sim_1_gt_p22_sim_1_imod_rec_resampled.txt')


#FSC_All_Now = np.loadtxt('FSC_All_gt_imod_resampled.txt')
#FSC_Meas_Now = np.loadtxt('FSC_In_gt_imod_resampled.txt')
#FSC_UnMeas_Now = np.loadtxt('FSC_Out_gt_imod_resampled.txt')


plt.plot(FSC_All_Now[:,0],FSC_All_Now[:,2])
plt.plot(FSC_Meas_Now[:,0],FSC_Meas_Now[:,2])
plt.plot(FSC_UnMeas_Now[:,0],FSC_UnMeas_Now[:,2])

plt.title('FSC IMod and GT')
plt.legend(['Usual FSC, IMOD vs GT','FSC in Measured Region','FSC in Not Measured Region'])
plt.xlim(0,0.2); plt.ylim(0,1.0)
plt.xlabel('Inverse Angstroms')
plt.ylabel('FSC')
plt.savefig('FSC_All_gt_imod_resampled.jpg')


#   Isonet
FSC_All_Now = np.loadtxt('FSC_All_p22_sim_1_gt_p22_sim_1_isonet_nomask_iter15.txt')
FSC_Meas_Now = np.loadtxt('FSC_In_p22_sim_1_gt_p22_sim_1_isonet_nomask_iter15.txt')
FSC_UnMeas_Now = np.loadtxt('FSC_Out_p22_sim_1_gt_p22_sim_1_isonet_nomask_iter15.txt')



plt.plot(FSC_All_Now[:,0],FSC_All_Now[:,2])
plt.plot(FSC_Meas_Now[:,0],FSC_Meas_Now[:,2])
plt.plot(FSC_UnMeas_Now[:,0],FSC_UnMeas_Now[:,2])

plt.title('FSC Isonet and GT')
plt.legend(['Usual FSC, Isonet vs GT','FSC in Measured Region','FSC in Not Measured Region'])
plt.xlim(0,0.2); plt.ylim(0,1.0)
plt.xlabel('Inverse Angstroms')
plt.ylabel('FSC')
plt.savefig('FSC_All_gt_isonet_nomask_not_resampled.jpg')




#   CN

FSC_All_Now = np.loadtxt('FSC_All_p22_sim_1_gt_p2s_1_cn_231128_lp18-c31_resampled.txt')
FSC_Meas_Now = np.loadtxt('FSC_In_p22_sim_1_gt_p2s_1_cn_231128_lp18-c31_resampled.txt')
FSC_UnMeas_Now = np.loadtxt('FSC_Out_p22_sim_1_gt_p2s_1_cn_231128_lp18-c31_resampled.txt')




plt.plot(FSC_All_Now[:,0],FSC_All_Now[:,2])
plt.plot(FSC_Meas_Now[:,0],FSC_Meas_Now[:,2])
plt.plot(FSC_UnMeas_Now[:,0],FSC_UnMeas_Now[:,2])

plt.title('FSC Crowther and GT')
plt.legend(['Usual FSC, Crowther vs GT','FSC in Measured Region','FSC in Not Measured Region'])
plt.xlim(0,0.2); plt.ylim(0,1.0)
plt.xlabel('Inverse Angstroms')
plt.ylabel('FSC')
plt.savefig('FSC_All_gt_p2s_1_cn_231128_lp18-c31_resampled.jpg')






#   eman


FSC_All_Now = np.loadtxt('FSC_All_p22_sim_1_gt_p22_sim_1_eman_new.txt')
FSC_Meas_Now = np.loadtxt('FSC_In_p22_sim_1_gt_p22_sim_1_eman_new.txt')
FSC_UnMeas_Now = np.loadtxt('FSC_Out_p22_sim_1_gt_p22_sim_1_eman_new.txt')




plt.plot(FSC_All_Now[:,0],FSC_All_Now[:,2])
plt.plot(FSC_Meas_Now[:,0],FSC_Meas_Now[:,2])
plt.plot(FSC_UnMeas_Now[:,0],FSC_UnMeas_Now[:,2])

plt.title('FSC EMAN and GT')
plt.legend(['Usual FSC, EMAN vs GT','FSC in Measured Region','FSC in Not Measured Region'])
plt.xlim(0,0.2); plt.ylim(0,1.0)
plt.xlabel('Inverse Angstroms')
plt.ylabel('FSC')
plt.savefig('FSC_All_gt_p22_sim_1_eman_new.jpg')







#%%      Section 4a Create FSC Outputs; n1 and n2 are the normalizations of
#  f and g. cH means the cosine of the phase residual.

nx2 = nx/2;
ny2 = ny/2;
nz2 = nz/2;

lsd2=nx+2;

dx2 = 1.0/float(nx2)/float(nx2);
dy2 = 1.0/float(ny2)/float(ny2);
dz2 = 1.0/float(nz2)/float(nz2);
# int inc = Util::round(float(std::max(std::max(nx2,ny2),nz2))/w);
w=1;

inc = max(nx2,ny2,nz2)/w;
inc = int(inc)




startTime = time.time()


#[retcHGlobal,lr]     = CreateFTLikeOutputs(inc,nx,ny,nz,dcH)
[ret,n1,n2,lr] = CreateFSCOutputs(inc,nx,ny,nz,d1,d2);
#[FPower,lr]    = CreateFTLikeOutputs(inc,nx,ny,nz,dFPower)


deltaTime =time.time()-startTime;
print("CreateFSCOutputs performed in %f seconds for size nx=%g " % (deltaTime,nx))

#%%

#%%      Section 4b  Write out FSCs. Define RMax based on this


# for nx =32, linc is 17
 
fNResRoot ='Test'
resultAveOut = fNResRoot+'globalFSC.csv';
APixels=1

linc = 0;
for i in range(inc+1):
    if (lr[i]>0):
        linc +=1;

result  = [0 for i in range(3*linc)];

ii = -1;
for i in range(inc+1):
    if (lr[i]>0): 
        ii +=1;
        result[ii]        = float(i)/float(2*inc);
        result[ii+linc]   = float(ret[i] / (np.sqrt(n1[i] * n2[i])));
        result[ii+2*linc] = lr[i]  ;# Number of Values


NormalizedFreq = result[0:(inc+1)];
resultAve= result[(inc+1):(2*(inc+1))];# This takes values inc+1 values from inc+1 to 2*inc+1

with open(resultAveOut, "w") as fL1:
    AveWriter = csv.writer(fL1)
    for j in range(inc+1):
        valFreqNormalized = NormalizedFreq[j];
        valFreq           = valFreqNormalized/APixels;
        valFSCshell= resultAve[j];
        AveWriter.writerow([valFreqNormalized,valFreq,valFSCshell])

#k0  = np.array(range(Nxf))/APixels/2.0/Nxf;


aa= np.abs(np.array(resultAve))<.13
bb=np.where(aa)[0];
try:
    RMax=bb[0]+4
except:
    RMax=inc
    
    
RMax= min(RMax,inc)    
print('simple FSC written out to '+resultAveOut)
print('RMax = %d'%(RMax))








#%%


if 1:
    h5f_HalfMap1= mrcfile.open(fNHalfMap1)
    f= h5f_HalfMap1.data.T

    h5f_HalfMap2= mrcfile.open(fNHalfMap2)
    g= h5f_HalfMap2.data.T

    h5f_HalfMap1.close()
    h5f_HalfMap2.close()


     
startTime = time.time()

#fT=f.T;# Now it is like EMAN
#gT=g.T;#



#[nx,ny,nz] =fT.shape
[nx,ny,nz] =f.shape

deltaTime =time.time()-startTime;
print("Maps read in %f seconds for size nx=%g " % (deltaTime,nx))





#%%      Section  Add Axes For Checking
#     
if 0:
    fTPlus=AddAxes(fT,2,10)
    
    h5f_write = h5py.File('fTPlus.hdf','w')
    h5f_write.create_dataset('MDF/images/0/image',data=fTPlus)
    # <HDF5 dataset "array": shape (63, 63, 63), type "<f8">
    h5f_write.close()



#%%      Section 2  Take Fourier Transform; find Cos phase residual
#                    CosHangle


startTime = time.time()


#[fp,gp]=ZeroPad(nx,ny,nz,fT,gT)
[fp,gp]=ZeroPad(nx,ny,nz,f,g)


deltaTime =time.time()-startTime;
print("NormPad created in %f seconds for size nx=%g " % (deltaTime,nx))


startTime = time.time()

F =np.fft.fftn(fp);#F= temp.transpose();
G =np.fft.fftn(gp);#G= temp.transpose();


H=F*np.conj(G);
H.shape
HAngle= np.angle(H);
CosHangle=np.cos(HAngle);


deltaTime =time.time()-startTime;
print("FFTs performed in %f seconds for size nx=%g " % (deltaTime,nx))

#FFTs performed in 4.559401 seconds for size nx=256 



#%%




# for nx =32, linc is 17


linc = 0;
for i in range(inc+1):
    if (lr[i]>0):
        linc +=1;

result  = [0 for i in range(3*linc)];

ii = -1;
for i in range(inc+1):
    if (lr[i]>0): 
        ii +=1;
        result[ii]        = float(i)/float(2*inc);
        result[ii+linc]   = float(ret[i] / (np.sqrt(n1[i] * n2[i])));
        result[ii+2*linc] = lr[i]  ;# Number of Values


NormalizedFreq = result[0:(inc+1)];
resultAve= result[(inc+1):(2*(inc+1))];# This takes values inc+1 values from inc+1 to 2*inc+1

with open(resultAveOut, "w") as fL1:
    AveWriter = csv.writer(fL1)
    for j in range(inc+1):
        valFreqNormalized = NormalizedFreq[j];
        valFreq           = valFreqNormalized/APixels;
        valFSCshell= resultAve[j];
        AveWriter.writerow([valFreqNormalized,valFreq,valFSCshell])

#k0  = np.array(range(Nxf))/APixels/2.0/Nxf;


aa= np.abs(np.array(resultAve))<.13
bb=np.where(aa)[0];
try:
    RMax=bb[0]+4
except:
    RMax=inc
    
    
RMax= min(RMax,inc)    
print('simple FSC written out to '+resultAveOut)
print('RMax = %d'%(RMax))












#%%

with open(resultAveOut, "w") as fL1:
    AveWriter = csv.writer(fL1)
    for j in range(inc+1):
        valFreqNormalized = NormalizedFreq[j];
        valFreq           = valFreqNormalized/APixels;
        valFSCshell= resultAve[j];
        AveWriter.writerow([valFreqNormalized,valFreq,valFSCshell])

#k0  = np.array(range(Nxf))/APixels/2.0/Nxf;


aa= np.abs(np.array(resultAve))<.13
bb=np.where(aa)[0];
try:
    RMax=bb[0]+4
except:
    RMax=inc
    
    
RMax= min(RMax,inc)    
print('simple FSC written out to '+resultAveOut)
print('RMax = %d'%(RMax))


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



def createFTarrays(nx,ny,nz,lsd2,lr,inc,dx2,dy2,dz2,dcH,dFPower):

    lrMaxOver2= int(lr[-1]//2);
               
    kXofR   = np.zeros([inc+1,lrMaxOver2],dtype=int)
    kYofR   = np.zeros([inc+1,lrMaxOver2],dtype=int)
    kZofR   = np.zeros([inc+1,lrMaxOver2],dtype=int)
    retcH   = np.zeros([inc+1,lrMaxOver2])
    retFT   = np.zeros([inc+1,lrMaxOver2])
    n12ofR  = np.zeros([inc+1,lrMaxOver2])
       
    NumAtEachR = np.zeros(inc+1,dtype=int);
#               
    rmax=0;
    for iz in range(nz):
        kz=iz;
        if (iz>nz2):
            kz=iz-nz;# This is the actual value kz, can be neg
        argz= float(kz*kz)*dz2;
        for iy in range(ny):
            ky=iy;# This is the actual value ky, can be neg
            if (iy>ny2):
                ky=iy-ny;
            argy= argz+float(ky*ky)*dy2;
            for ix in range(0,nx,2):
                if ( (ix==0) & (ky<0)): continue
                if ( (ix==0) & (ky==0)& (kz<0)): continue
                kx=float(ix)/2.0;
                argx = 0.5*np.sqrt(argy + kx*kx*dx2);
                r=int(round(inc*2*argx))
                if(r <= inc):
                    NumAtEachR[r]+=1;
                    LastInd = NumAtEachR[r]-1;
                    #ii = ix + (iy  + iz * ny)* lsd2;
                    kXofR[r,LastInd] =int(round(kx));
                    kYofR[r,LastInd] =ky;
                    kZofR[r,LastInd] =kz;
                    retcHNow  = dcH[ix  ,iy,iz];
                    retFTNow  = dFPower[ix  ,iy,iz];
                    n12rNow   = 1;
                     #retofRR[r]=retNowR;  retofRI[r]=retNowI;  
                     #n1ofR[r]=n1Now;   n2ofR[r]=n2Now;
                    retcH[r, LastInd]=retcHNow;
                    retFT[r, LastInd]=retFTNow;
                    n12ofR[r, LastInd]=n12rNow;   

                    if r>rmax: rmax=r;
    print(rmax)
    return [kXofR, kYofR, kZofR,retcH,retFT,n12ofR]

    
    
def createFSCarrays(nx,ny,nz,lsd2,lr,inc,dx2,dy2,dz2,d1,d2):

    lrMaxOver2= int(lr[-1]//2);
               
    kXofR   = np.zeros([inc+1,lrMaxOver2],dtype=int)
    kYofR   = np.zeros([inc+1,lrMaxOver2],dtype=int)
    kZofR   = np.zeros([inc+1,lrMaxOver2],dtype=int)
    retofRR = np.zeros([inc+1,lrMaxOver2])
    retofRI = np.zeros([inc+1,lrMaxOver2])
    n1ofR   = np.zeros([inc+1,lrMaxOver2])
    n2ofR   = np.zeros([inc+1,lrMaxOver2])
       
    NumAtEachR = np.zeros(inc+1,dtype=int);
#               
    rmax=0;
    for iz in range(nz):
        kz=iz;
        if (iz>nz2):
            kz=iz-nz;# This is the actual value kz, can be neg
        argz= float(kz*kz)*dz2;
        for iy in range(ny):
            ky=iy;# This is the actual value ky, can be neg
            if (iy>ny2):
                ky=iy-ny;
            argy= argz+float(ky*ky)*dy2;
            for ix in range(0,nx,2):
                if ( (ix==0) & (ky<0)): continue
                if ( (ix==0) & (ky==0)& (kz<0)): continue
                kx=float(ix)/2.0;
                argx = 0.5*np.sqrt(argy + kx*kx*dx2);
                r=int(round(inc*2*argx))
                if(r <= inc):
                    NumAtEachR[r]+=1;
                    LastInd = NumAtEachR[r]-1;
                    #ii = ix + (iy  + iz * ny)* lsd2;
                    kXofR[r,LastInd] =int(round(kx));
                    kYofR[r,LastInd] =ky;
                    kZofR[r,LastInd] =kz;
                    retrRNow  = d1[ix  ,iy,iz] * d2[ix  ,iy,iz];
                    retrRNow += d1[ix+1,iy,iz] * d2[ix+1,iy,iz];
                    retrINow  = d1[ix  ,iy,iz] * d2[ix+1,iy,iz];
                    retrINow -= d1[ix+1,iy,iz] * d2[ix  ,iy,iz];
                    n1rNow    = d1[ix  ,iy,iz] * d1[ix  ,iy,iz];
                    n1rNow   += d1[ix+1,iy,iz] * d1[ix+1,iy,iz];
                    n2rNow    = d2[ix  ,iy,iz] * d2[ix  ,iy,iz];
                    n2rNow   += d2[ix+1,iy,iz] * d2[ix+1,iy,iz];
                     #retofRR[r]=retNowR;  retofRI[r]=retNowI;  
                     #n1ofR[r]=n1Now;   n2ofR[r]=n2Now;
                    retofRR[r, LastInd]=retrRNow;
                    retofRI[r, LastInd]=retrINow;
                    n1ofR[r, LastInd]=n1rNow;   
                    n2ofR[r, LastInd]=n2rNow;

                    if r>rmax: rmax=r;
    print(rmax)
    return [kXofR, kYofR, kZofR,retofRR,retofRI,n1ofR,n2ofR, NumAtEachR]
