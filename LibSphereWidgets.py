#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 16:16:52 2019

@author: pbaldwin
"""

import numpy as np
from sys import argv
import csv
import sys, os
import pandas as pd
import matplotlib
import scipy
from scipy import special
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import glob;
import time
import shutil;
#from math import *
#from numba import *
#  from numbapro import cuda


#%%

from numba import autojit
#@autojit    
def AnglesToVecsTetrahedral(AngleTiltsOrNVecs0,AngleRotsOrFlag,CylinderLengthVec=0):#DSym=1 or 2

    # The  Relion asymmetric unit is a f f' on Fig 2 of the Transform paper
    # where f' is the centroid of a neighboring cell

    # CylinderLengthVec is the length of the Cylinder if a chimera bild file is used
    
    lvl0=0;         # There is a face along z
    lvl1=109.4712;  #  that is acos(-1/3)  ; There  are 3 faces at this angle


    TET = np.array([0,lvl0,0,   0,lvl0,120,    0,lvl0,240,
      0,lvl1,60,   0,lvl1,180,    0,lvl1,300,
      120,lvl1,60, 120,lvl1,180,  120,lvl1,300,
      240,lvl1,60, 240,lvl1,180,  240,lvl1,300 ])

    Skip = 12; # for Octahedron

    if 0:
        TET = np.array([0,lvl0,0,  
                  0,lvl1,60,   0,lvl1,180,    0,lvl1,300 ])

        Skip = 4; # for Octahedron

##


    AllTETMat= np.zeros([Skip,3,3])
    TETtriples = np.reshape(TET,(Skip,3))
    
    for jTET,TETNow in enumerate(TETtriples):
        # EMAN variables
        azITET  = TETNow[0]*np.pi/180.0;
        altTET = TETNow[1]*np.pi/180.0;
        phiTET = TETNow[2]*np.pi/180.0;
        Zaz  = np.matrix([[np.cos(azITET), np.sin(azITET), 0],  [-np.sin(azITET), np.cos(azITET), 0],  [0, 0, 1]])
        Xalt = np.matrix([ [1, 0, 0], [0,np.cos(altTET), np.sin(altTET)],  [0,-np.sin(altTET), np.cos(altTET)] ])
        Zphi = np.matrix([[np.cos(phiTET), np.sin(phiTET), 0],  [-np.sin(phiTET), np.cos(phiTET), 0],  [0, 0, 1]])
        #AllThreeNow = Zaz*Xalt*Zphi
        AllThreeNow = Zphi*Xalt*Zaz; # This is B version
        #AllThreeNow = Xalt*Zphi
        #AllThreeNow = Zaz*Xalt; # Bad
        AllTETMat[jTET,:,:] = AllThreeNow;


    if 0:
        for jAngle in range(sAngleTilts):
            #jj = np.random.randint(1,8000+1)
            jj = jAngle
            AngleRotNow  = AngleRots[jj];
            AngleTiltNow = AngleTilts[jj];
            #LengthNow= CylinderLengthVec[jj]
            c1= np.cos(AngleRotNow)
            s1= np.sin(AngleRotNow)
            cy= np.cos(AngleTiltNow)
            sy= np.sin(AngleTiltNow)
            nx= c1*sy;
            ny= s1*sy;
            nz= cy;
            NVecs0[:,jAngle] = [nx,ny,nz]
    
    if len(np.shape(AngleRotsOrFlag))==0:
        NVecs0= AngleTiltsOrNVecs0.T
        sAngleTilts=np.shape(NVecs0)[1]
    else:
        AngleTilts = AngleTiltsOrNVecs0
        sAngleTilts=len(AngleTilts);
        AngleRots = AngleRotsOrFlag
        nx= np.cos(AngleRots)*np.sin(AngleTilts)
        ny= np.sin(AngleRots)*np.sin(AngleTilts)
        nz=                   np.cos(AngleTilts)
        NVecs0= np.hstack((nx,ny,nz)).reshape(3,sAngleTilts)
        
    #sAngleTilts= 4;
    #NVecs0          = np.zeros([3,sAngleTilts]);
    NVecs           = np.matrix(np.zeros([sAngleTilts*Skip,3]));
    CylinderVecSymm = np.zeros(sAngleTilts*Skip);


    zhat = np.array([0,0,1]).T
    TestVertices =  np.dot(AllTETMat, zhat)
    NVecs0AndCopies = np.dot(AllTETMat, NVecs0)

    for jAngle in range(sAngleTilts):
        for jCopy in range(Skip):   
            NVecs[jAngle*Skip+jCopy,:] = NVecs0AndCopies[jCopy,:,jAngle]


    return NVecs


#%%    Section 0;  routine 2;  AnglesToVecsCDOctahedral for  AnglesToVecs  

#  AngleTilts = ThetaVec;
#  AngleRots =PhiVec;

#  AngleTilts = ThetaVec;
#  AngleRots =PhiVec;

from numba import autojit
#@autojit    
def AnglesToVecsOctahedral(AngleTiltsOrNVecs0,AngleRotsOrFlag,CylinderLengthVec=0):#DSym=1 or 2

    # The  Relion asymmetric unit is a f f' on Fig 2 of the Transform paper
    # where f' is the centroid of a neighboring cell

    # CylinderLengthVec is the length of the Cylinder if a chimera bild file is used
    lvl0=0.;
    lvl1=90.;
    lvl2=180.;
    OCT = np.array([0,lvl0,0,   0,lvl0,90,    0,lvl0,180,    0,lvl0,270,
      0,lvl1,0,   0,lvl1,90,    0,lvl1,180,    0,lvl1,270,
      90,lvl1,0,  90,lvl1,90,   90,lvl1,180,   90,lvl1,270,
      180,lvl1,0, 180,lvl1,90,  180,lvl1,180,  180,lvl1,270,
      270,lvl1,0, 270,lvl1,90,  270,lvl1,180,  270,lvl1,270,
      0,lvl2,0,   0,lvl2,90,    0,lvl2,180,    0,lvl2,270 ])
    Skip = 24; # for Octahedron

    if 0:
        OCT = np.array([0,lvl0,0,   0,lvl0,90,    0,lvl0,180,    0,lvl0,270,
          0,lvl1,0,   0,lvl1,90,    0,lvl1,180,    0,lvl1,270,
          90,lvl1,0,  90,lvl1,90,   90,lvl1,180,   90,lvl1,270 ])
        Skip = 12; # for Octahedron
    
    
    if 0:   
        AdjustTemp=0
        OCT = np.array([
              0,lvl0,AdjustTemp,                           # z -> z
              0,lvl1,AdjustTemp,   0,lvl1,AdjustTemp+90,    0,lvl1,AdjustTemp+180,    0,lvl1,AdjustTemp+270 , # z ->+- x; z -> +- y
              0,lvl2,AdjustTemp])      # # z -> -z
    
        Skip = 6; # for Octahedron


    AllOCTMat= np.zeros([Skip,3,3])
    OCTtriples = np.reshape(OCT,(Skip,3))
    
    for jOCT,OCTNow in enumerate(OCTtriples):
        # EMAN variables
        azIOCT  = OCTNow[0]*np.pi/180.0;
        altOCT = OCTNow[1]*np.pi/180.0;
        phiOCT = OCTNow[2]*np.pi/180.0;
        Zaz  = np.matrix([[np.cos(azIOCT), np.sin(azIOCT), 0],  [-np.sin(azIOCT), np.cos(azIOCT), 0],  [0, 0, 1]])
        Xalt = np.matrix([ [1, 0, 0], [0,np.cos(altOCT), np.sin(altOCT)],  [0,-np.sin(altOCT), np.cos(altOCT)] ])
        Zphi = np.matrix([[np.cos(phiOCT), np.sin(phiOCT), 0],  [-np.sin(phiOCT), np.cos(phiOCT), 0],  [0, 0, 1]])
        #AllThreeNow = Zaz*Xalt*Zphi
        AllThreeNow = Zphi*Xalt*Zaz; # This is B version
        #AllThreeNow = Xalt*Zphi
        #AllThreeNow = Zaz*Xalt; # Bad
        AllOCTMat[jOCT,:,:] = AllThreeNow;

    
    if len(np.shape(AngleRotsOrFlag))==0:
        NVecs0= AngleTiltsOrNVecs0.T
        sAngleTilts=np.shape(NVecs0)[1]
    else:
        AngleTilts = AngleTiltsOrNVecs0
        sAngleTilts=len(AngleTilts);
        AngleRots = AngleRotsOrFlag
        nx= np.cos(AngleRots)*np.sin(AngleTilts)
        ny= np.sin(AngleRots)*np.sin(AngleTilts)
        nz=                   np.cos(AngleTilts)
        NVecs0= np.hstack((nx,ny,nz)).reshape(3,sAngleTilts)
        
    #sAngleTilts= 4;
    #NVecs0          = np.zeros([3,sAngleTilts]);
    NVecs           = np.matrix(np.zeros([sAngleTilts*Skip,3]));
    CylinderVecSymm = np.zeros(sAngleTilts*Skip);

    NVecs0AndCopies = np.dot(AllOCTMat, NVecs0)

    for jAngle in range(sAngleTilts):
        for jCopy in range(Skip):   
            NVecs[jAngle*Skip+jCopy,:] = NVecs0AndCopies[jCopy,:,jAngle]


    return NVecs



#%%    Section 0;  routine 2;  AnglesToVecsIcosahedral for  AnglesToVecs  

#  AngleTilts = ThetaVec;
#  AngleRots =PhiVec;

#  AngleTilts = ThetaVec;
#  AngleRots =PhiVec;

from numba import autojit
#@autojit    
def AnglesToVecsIcosahedral(AngleTiltsOrNVecs0,AngleRotsOrFlag,CylinderLengthVec=0):#
#def AnglesToVecsIcosahedral(AngleTilts,AngleRots):#DSym=1 or 2

    # The  Relion asymmetric unit is a f f' on Fig 2 of the Transform paper
    # where f' is the centroid of a neighboring cell

    # CylinderLengthVec is the length of the Cylinder if a chimera bild file is used
    lvl0 = 0.0; # there is one pentagon on top; five-fold along z
    lvl1 = (180/np.pi)*np.arctan(2.0);# 63.4349; // that is atan(2)  // there are 5 pentagons with centers at this height (angle)
    lvl2 = 180.0 - lvl1;#116.5651; //that is 180-lvl1  // there are 5 pentagons with centers at this height (angle)
    lvl3 = 180.0;
    
    ICOS = np.array([  # This is a pentagon  normal to z
              0,lvl0,0,    0,lvl0,288,   0,lvl0,216,   0,lvl0,144,  0,lvl0,72,
              0,lvl1,36,   0,lvl1,324,   0,lvl1,252,   0,lvl1,180,  0,lvl1,108,
              72,lvl1,36,  72,lvl1,324,  72,lvl1,252,  72,lvl1,180,  72,lvl1,108,
              144,lvl1,36, 144,lvl1,324, 144,lvl1,252, 144,lvl1,180, 144,lvl1,108,
              216,lvl1,36, 216,lvl1,324, 216,lvl1,252, 216,lvl1,180, 216,lvl1,108,
              288,lvl1,36, 288,lvl1,324, 288,lvl1,252, 288,lvl1,180, 288,lvl1,108,
               36,lvl2,0,   36,lvl2,288,  36,lvl2,216,  36,lvl2,144,  36,lvl2,72,
              108,lvl2,0,  108,lvl2,288, 108,lvl2,216, 108,lvl2,144, 108,lvl2,72,
              180,lvl2,0,  180,lvl2,288, 180,lvl2,216, 180,lvl2,144, 180,lvl2,72,
              252,lvl2,0,  252,lvl2,288, 252,lvl2,216, 252,lvl2,144, 252,lvl2,72,
              324,lvl2,0,  324,lvl2,288, 324,lvl2,216, 324,lvl2,144, 324,lvl2,72,
                0,lvl3,0,    0,lvl3,288,   0,lvl3,216,   0,lvl3,144,   0,lvl3,72])
   
    Skip = 60; # for Octahedron

    if 0:
        ICOS = np.array([  # This is a pentagon  normal to z
              0,lvl0,0,  
              0,lvl1,36,   0,lvl1,324,   0,lvl1,252,   0,lvl1,180,  0,lvl1,108,
               36,lvl2,0,   36,lvl2,288,  36,lvl2,216,  36,lvl2,144,  36,lvl2,72,
                0,lvl3,0])

        Skip = 12; # for Octahedron

    AllICOSMat= np.zeros([Skip,3,3])
    ICOStriples = np.reshape(ICOS,(Skip,3))
    
    for jICOS,ICOSNow in enumerate(ICOStriples):
        # EMAN variables
        azICOS  = ICOSNow[0]*np.pi/180.0;
        altICOS = ICOSNow[1]*np.pi/180.0;
        phiICOS = ICOSNow[2]*np.pi/180.0;
        Zaz  = np.matrix([[np.cos(azICOS), np.sin(azICOS), 0],  [-np.sin(azICOS), np.cos(azICOS), 0],  [0, 0, 1]])
        Xalt = np.matrix([ [1, 0, 0], [0,np.cos(altICOS), np.sin(altICOS)],  [0,-np.sin(altICOS), np.cos(altICOS)] ])
        Zphi = np.matrix([[np.cos(phiICOS), np.sin(phiICOS), 0],  [-np.sin(phiICOS), np.cos(phiICOS), 0],  [0, 0, 1]])
        AllThreeNow = Zaz*Xalt*Zphi
        AllThreeNow = Zphi*Xalt*Zaz; # This is B version
        #AllThreeNow = Xalt*Zphi
        #AllThreeNow = Zaz*Xalt; # Bad
        AllICOSMat[jICOS,:,:] = AllThreeNow;

    #print(np.shape(AllICOSMat))


     
    if len(np.shape(AngleRotsOrFlag))==0:
        NVecs0= AngleTiltsOrNVecs0.T
        sAngleTilts=np.shape(NVecs0)[1]
    else:
        AngleTilts = AngleTiltsOrNVecs0
        sAngleTilts=len(AngleTilts);
        AngleRots = AngleRotsOrFlag
        nx= np.cos(AngleRots)*np.sin(AngleTilts)
        ny= np.sin(AngleRots)*np.sin(AngleTilts)
        nz=                   np.cos(AngleTilts)
        NVecs0= np.hstack((nx,ny,nz)).reshape(3,sAngleTilts)
        
    #sAngleTilts= 4;
    #NVecs0          = np.zeros([3,sAngleTilts]);
    NVecs           = np.matrix(np.zeros([sAngleTilts*Skip,3]));
    CylinderVecSymm = np.zeros(sAngleTilts*Skip);

  
    NVecs0AndCopies = np.dot(AllICOSMat, NVecs0)

    for jAngle in range(sAngleTilts):
        for jCopy in range(Skip):   
            NVecs[jAngle*Skip+jCopy,:] = NVecs0AndCopies[jCopy,:,jAngle]

        
#        for id in range(Skip):
#            Z1AngleNow = ICOS[id*3]   * np.pi/180.0;
#            XAngleNow  = ICOS[id*3+1] * np.pi/180.0;
#            Z2AngleNow = ICOS[id*3+2] * np.pi/180.0;
#            c= np.cos(Z1AngleNow);s=np.sin(Z1AngleNow)
#            Z1Now= np.matrix([[c,s,0],[-s,c,0],[0,0,1]])
#            c= np.cos(XAngleNow);s=np.sin(XAngleNow)
#            XNow= np.matrix([[1,0,0],[0,c,s],[0,-s,c]])
#            c= np.cos(Z2AngleNow);s=np.sin(Z2AngleNow)
#            Z2Now= np.matrix([[c,s,0],[-s,c,0],[0,0,1]])
#            RotNow = Z1Now*XNow*Z2Now;
#            IndexOut = jAngle*Skip+id;
#            NVecs[IndexOut,:]=NVec0*RotNow
#            #CylinderVecSymm[IndexOut]=LengthNow;

    #return NVecs,CylinderVecSymm
    return NVecs 

  

#%%     Section 0;   routine 3; Convert Angles to Vecs;  with C or D symmetries
        
from numba import autojit
#@autojit    
def AnglesToVecsCD(AngleTiltsOrNVecs0,AngleRotsOrFlag,CSym,CorD):#CorD=1 or 2


    if len(np.shape(AngleRotsOrFlag))==0:
        NVecs0T= AngleTiltsOrNVecs0
        AngleTilts =np.arccos(NVecs0T[:,2])*180/np.pi
        AngleRots = np.arctan2(NVecs0T[:,1], NVecs0T[:,0])*180/np.pi
        #print(AngleTilts,AngleRots)
    else:
        AngleTilts = AngleTiltsOrNVecs0
        AngleRots =  AngleRotsOrFlag;
        
    sAngleTilts=len(AngleTilts);# am not sure what is the difference with AnglesToVecsCDv2
    Skip= CorD*CSym;
    NVecs = np.zeros([sAngleTilts*Skip,3])
    for jAngle in range(sAngleTilts):
        AngleRotNow0  = AngleRots[jAngle];
        AngleTiltNow = AngleTilts[jAngle];
        for jAngSym in range(CSym):
            jAngAdd = jAngSym*360.0/CSym;
            AngleRotNow  = AngleRotNow0+jAngAdd;
            for jDSym in range(CorD):
                if jDSym==1: AngleTiltNow = 180.0- AngleTiltNow;
                c1= np.cos(AngleRotNow*np.pi/180.0)
                s1= np.sin(AngleRotNow*np.pi/180.0)
                cy= np.cos(AngleTiltNow*np.pi/180.0)
                sy= np.sin(AngleTiltNow*np.pi/180.0)
                nx= c1*sy;
                ny= s1*sy;
                nz= cy;
                if np.isnan(nx): continue
                NVecs[jAngle*Skip + jAngSym*CorD + jDSym,:]=[nx,ny,nz]       
    return NVecs



#%%

def AngleFile2nVec(AngleFileName):
    OutPutAngles = np.loadtxt(AngleFileName)
    NumAngles = np.shape(OutPutAngles)[0]
    #ZerosVec = np.zeros(NumAngles)
    psiVec   = OutPutAngles[:,0]*np.pi/180.0;
    thetaVec = OutPutAngles[:,1]*np.pi/180.0;
    phiVec  = OutPutAngles[:,2]*np.pi/180.0;

    nVecx = np.sin(thetaVec)*np.cos(phiVec)
    nVecy = np.sin(thetaVec)*np.sin(phiVec)
    nVecz = np.cos(thetaVec)
    #The angle file, FileNameOfAngles, should be formatted as three columns for the Z Y Z convention:
    #that is 'PSI', 'THETA', 'PHI'.
    #The first column is irrelevant for the purposes of the sampling. It is PSI. The second column is the tilts: which probably is positive and ranges from 0 to 180.

    nVec = np.vstack((nVecx, nVecy,nVecz)).T
    return nVec,OutPutAngles
# nVec =  AngleFile2nVec(AngleFileName)
    
    return nVec

#%%
#@autojit
def GetIFibDirections_Now(forceNumberofRaystothisinteger):
    #deltaTheta = float(argv[1]); 
    #deltaThetaRad=np.pi*deltaTheta/180.0;
    #Directions_csv_out = argv[2];
    
    # Now, deltaTheta takes up a cone which has area
    #    2 pi * (1-cos(deltaTheta))
    # The surface area of the Hemisphere is 2 pi
    # We need Hemisphere (as opposed to sphere) because the 
    #   unit vectors take up both directions
    #  So the fraction for each plaquette is (1- cos(theta))
    # N should be about 1 over this fraction.
    # But we will create twice so many.
    
    
    #  Section 1. Calculate location of points
    #
    #FractionOfHemi = 1- np.cos(deltaThetaRad);
    #NumPoints = int(round(2.0/FractionOfHemi));
    NumPoints = forceNumberofRaystothisinteger;
    
    phi=(1.0+ np.sqrt(5.0))/2.0-1
    ga = phi*2.0*np.pi;#  The Golden Angle
    
    
    nVec=np.zeros([NumPoints,3]);
    jLongitude=np.arange(NumPoints);
    
    longitude = ga*jLongitude;
    latitude  =  np.arcsin(1.0*jLongitude/NumPoints);# for hemisphere
    nhatz = np.sin(latitude);
    nhatx = np.cos(latitude)*np.cos(longitude);
    nhaty = np.cos(latitude)*np.sin(longitude);

    nVec = np.concatenate([nhatx,nhaty,nhatz])
    nVec = np.reshape(nVec,[3,NumPoints])

    return nVec.T


#%%
#@autojit
def ReturnSampling(SpiralVec,AllProjPoints,FourierRadius,TiltInDeg,NumberForEachTilt):
    
    tiltRad = TiltInDeg*np.pi/180;
    
    # AllProjPoints may be very large, so we need to chunk it
    InnerProductDirMatrix=np.dot(SpiralVec ,AllProjPoints.T);

    vv=np.where(np.abs(InnerProductDirMatrix)*FourierRadius<0.5)
    #print(type(vv))

    InnerProdBoolean = np.zeros_like(InnerProductDirMatrix)
    InnerProdBoolean[vv[0],vv[1]]=1

    InnerProdBooleanSum0 = np.sum(InnerProdBoolean,axis=0)
    InnerProdBooleanSum1 = np.sum(InnerProdBoolean,axis=1)

    if TiltInDeg==0:    
        return InnerProdBooleanSum0,InnerProdBooleanSum1,AllProjPoints
    
    print('Tilted')
    
    AllPointsPerp=  np.argmin(np.abs(AllProjPoints),axis=1)
    print('np.shape(AllProjPoints), np.shape(AllPointsPerp)')
    print(np.shape(AllProjPoints), np.shape(AllPointsPerp))
    #plt.hist(AllPointsPerp)

    PerpVector = np.zeros_like(AllProjPoints)
    PerpPerpVector = np.zeros_like(AllProjPoints)

    # Perp0
    Perp0 = np.where(AllPointsPerp==0)[0]

    nx =  AllProjPoints[Perp0,0]; nx =np.asarray(nx)
    ny =  AllProjPoints[Perp0,1]; ny =np.asarray(ny)
    nz =  AllProjPoints[Perp0,2]; nz =np.asarray(nz)
    
    normyz = np.sqrt(ny*ny+nz*nz)
    
    PerpVector[Perp0,1]= -nz/normyz
    PerpVector[Perp0,2]=  ny/normyz


    PerpPerpVector[Perp0,0]= nz*nz/normyz; 
    PerpPerpVector[Perp0,0]+=ny*ny/normyz
    PerpPerpVector[Perp0,1]=-nx*ny/normyz
    PerpPerpVector[Perp0,2]=-nx*nz/normyz; 

    # Perp1

    Perp1 = np.where(AllPointsPerp==1)[0]
    
    nx =  AllProjPoints[Perp1,0]; nx =np.asarray(nx)
    ny =  AllProjPoints[Perp1,1]; ny =np.asarray(ny)
    nz =  AllProjPoints[Perp1,2]; nz =np.asarray(nz)
    
    normxz = np.sqrt(nx*nx+nz*nz)
    PerpVector[Perp1,2]= -nx/normxz
    PerpVector[Perp1,0]=  nz/normxz

    PerpPerpVector[Perp1,0]=-nx*ny/normxz; 
    PerpPerpVector[Perp1,1]= nx*nx/normxz
    PerpPerpVector[Perp1,1]+=nz*nz/normxz
    PerpPerpVector[Perp1,2]=-ny*nz/normxz; 

    # Perp2

    Perp2 = np.where(AllPointsPerp==2)[0]
    
    nx =  AllProjPoints[Perp2,0]; nx =np.asarray(nx)
    ny =  AllProjPoints[Perp2,1]; ny =np.asarray(ny)
    nz =  AllProjPoints[Perp2,2]; nz =np.asarray(nz)
    
    normxy = np.sqrt(nx*nx+ny*ny)
    
    PerpVector[Perp2,0]= -ny/normxy
    PerpVector[Perp2,1]=  nx/normxy

    #PerpVector_N = preprocessing.normalize(PerpVector, norm='l2')

    PerpPerpVector[Perp2,0]=-nx*nz/normxy; 
    PerpPerpVector[Perp2,1]=-ny*nz/normxy
    PerpPerpVector[Perp2,2]= nx*nx/normxy
    PerpPerpVector[Perp2,2]+=ny*ny/normxy; 

    #PerpPerpVector_N = preprocessing.normalize(PerpPerpVector, norm='l2')
          
    NumPhi = int(90.*np.sin(tiltRad))
    NumAll = len(AllProjPoints)
    phi = np.linspace(0,2*np.pi,NumPhi)
    AllProjPointsAlong = np.cos(tiltRad) * AllProjPoints;
    AllProjPointsPerp  = np.sin(tiltRad) * PerpVector;
    AllProjPointsPP    = np.sin(tiltRad) * PerpPerpVector;


    AllProjAndTilts = np.zeros([NumAll*NumPhi,3])
    print(len(AllProjPoints))

    for jAll in range(NumAll):
        PerpNow = AllProjPointsPerp[jAll]
        PPNow   = AllProjPointsPP[jAll]
        PerpAllNow = np.outer(np.cos(phi),PerpNow)+ np.outer(np.sin(phi),PPNow)
        AllProjAndTilts[(jAll*NumPhi):((jAll+1)*NumPhi),:] = AllProjPointsAlong[jAll,:] +PerpAllNow  

    InnerProductDirMatrix=np.dot(SpiralVec ,AllProjAndTilts.T);

    vv=np.where(np.abs(InnerProductDirMatrix)*FourierRadius<0.5)
    #print(type(vv))

    InnerProdBoolean = np.zeros_like(InnerProductDirMatrix)
    InnerProdBoolean[vv[0],vv[1]]=1

    InnerProdBooleanSum0 = np.sum(InnerProdBoolean,axis=0)
    InnerProdBooleanSum1 = np.sum(InnerProdBoolean,axis=1)

    
    #ww= np.where(InnerProdBooleanSum1==0)[0]
    #NumZerosVec[jFR]  = len(ww)
    #MeanValueVec[jFR] = np.mean(InnerProdBooleanSum1)

    #print(FourierRadius,2*np.pi*FourierRadius*FourierRadius, NumFib)
    
    return InnerProdBooleanSum0,InnerProdBooleanSum1,AllProjAndTilts	

#%%
    
def ReturnSamplingLoop(SpiralVec,AllProjPoints,FourierRadius,TiltInDeg,NumberForEachTilt):
    
    sSpiralVec= len(SpiralVec);
    InnerProdBooleanSum1Sum= np.zeros(sSpiralVec)
    InnerProdBooleanSum0Sum= np.shape(AllProjPoints)[0]
    sAllProjPoints = len(AllProjPoints);
    
    ChunkSize = np.int(10000*3000/sSpiralVec);
    
    NumChunks= sAllProjPoints//ChunkSize
    print('sAllProjPoints = %g; '%(sAllProjPoints))
    print('NumChunks = %g; ChunkSize=%g'%(NumChunks,ChunkSize))
    
    #if NumChunks==0:
    jChunk=0
        
    for jChunk in range(NumChunks):
        AllProjPointsInner=AllProjPoints[jChunk*ChunkSize:(jChunk+1)*ChunkSize]
        InnerProdBooleanSum0,InnerProdBooleanSum1,AllProjAndTilts = \
            ReturnSampling(SpiralVec,AllProjPointsInner,FourierRadius,TiltInDeg,NumberForEachTilt)
        InnerProdBooleanSum1Sum += InnerProdBooleanSum1
        #InnerProdBooleanSum0Sum[jChunk*ChunkSize:(jChunk+1)*ChunkSize] = InnerProdBooleanSum0
        
    AllProjPointsInner = AllProjPoints[jChunk*ChunkSize:]
    InnerProdBooleanSum0,InnerProdBooleanSum1,AllProjAndTilts = \
         ReturnSampling(SpiralVec,AllProjPointsInner,FourierRadius,TiltInDeg,NumberForEachTilt)
    InnerProdBooleanSum1Sum += InnerProdBooleanSum1
    #InnerProdBooleanSum0Sum[jChunk*ChunkSize:] = InnerProdBooleanSum0

    return InnerProdBooleanSum0Sum,InnerProdBooleanSum1Sum,AllProjAndTilts	
