#!/usr/bin/env python

import os,sys,glob,time,copy
from math import floor,ceil,log,exp
import numpy as np
import dicom
from ROOT import TFile,TH3S,TH3D
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import arange, array, exp

def extrap1d(interpolator): 
    xs = interpolator.x
    ys = interpolator.y
    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)
    def ufunclike(xs):
      return array(map(pointwise, array(xs)))

    return ufunclike


# ==================================================================
# Configuration
# ==================================================================
A = {}
A['DATA_PATH'] = sys.argv[1]
A['PWD'] = os.environ['PWD']
# ==================================================================
# Get the plan info
# ==================================================================

print A['DATA_PATH'],A['PWD']
print A['PWD']
os.chdir(A['DATA_PATH'])
try: os.mkdir('out')
except:pass
os.chdir('out')
# ==================================================================
# Phantom
# ==================================================================
# Get CT data
#voxel_size = 8
f = TFile('phantom.root','recreate')
fCT = glob.glob(A['PWD']+'/'+A['DATA_PATH']+'/*.dcm')#change to .ima /.dcm
fCT.sort()
CT = dicom.read_file(fCT[-1])

format_str = '%sint%d' % (('u', '')[CT.PixelRepresentation],CT.BitsAllocated)


N = [0,0,0]
N[0] = CT.Rows
N[1] = CT.Columns
N[2] = len(fCT)
#  p = np.fromstring(CT.PixelData,np.dtype(format_str)).reshape((N[1],N[0]))

if not CT.SliceThickness==0: pix = CT.PixelSpacing+[CT.SliceThickness]
else: 
  try :pix = CT.PixelSpacing + [abs(dicom.read_file(fCT[0]).SliceLocation-dicom.read_file(fCT[1]).SliceLocation)]
  except :pix = CT.PixelSpacing + [abs(dicom.read_file(fCT[0]).ImagePositionPatient[2]-dicom.read_file(fCT[1]).ImagePositionPatient[2])]
    

pos = CT.ImagePositionPatient
#print pos
# Find the minimum z position

pos[2] = min( [dicom.read_file(dcm).ImagePositionPatient[2] for dcm in fCT] )
max = max( [dicom.read_file(dcm).ImagePositionPatient[2] for dcm in fCT] )
min = pos[2]
pix[2] = (max-min)/(len(fCT)-1)

pix=map(float,pix)
pos=map(float,pos)
#print pix,pos

# Clip edges against CT scan
h   = TH3S('hu','',N[0],pos[0],pos[0]+(N[0]-1)*pix[0],N[1],pos[1],pos[1]+(N[1]-1)*pix[1],N[2],pos[2],pos[2]+(N[2]-1)*pix[2])

rho = TH3D('rho','',N[0],pos[0],pos[0]+(N[0]-1)*pix[0],N[1],pos[1],pos[1]+(N[1]-1)*pix[1],N[2],pos[2],pos[2]+(N[2]-1)*pix[2])
# Cropped phantom
f.cd()
#print 'Voxel size (mm): original(%.3f,%.3f,%.3f), phantom(%.3f,%.3f,%.3f)' % tuple(pix+[voxel_size,voxel_size,voxel_size])

a = np.zeros((N[2],N[1],N[0]),format_str)
#HUUnit  = np.array([-4000,-1000, -727.9,-522.6       , 0            , 225.80    , 455.8 , 798   , 1203.1, 2500,3000,4000])  # Hounsfield Unit
HUUnit = np.array([-1027,-1000,-533.88,-99.13,-44.50,3.50,38.75,73.88,240.38,250.38,481.38,881.50,1310.13,2000.0, 2500.0, 3071.0])
RhoUnit = np.array([0.00324,0.00324,0.444,0.925,0.965,1.0,1.02,1.058,1.099,1.105,1.278,1.47,1.69,2.026,2.28,2.575]) #electron density
## Vacuum air LN450 Breast water Muscle Brain InnerBone Bone200 CB230% CB250% CorticalBone 
#RhoUnit = np.array([0.000001,0.0001, 0.292  , 0.438  , 1.000        , 1.099     , 1.285 , 1.473 , 1.707,  2.66497655,2.8419,2.8419])  #Electron density
rhoFit  = interp1d(HUUnit,RhoUnit)

z = 0
for dcm in fCT:
  CT = dicom.read_file(dcm)
  #print CT
  #print CT.InstanceNumber
  #z=CT.InstanceNumber-1
  #print z, CT.ImagePositionPatient[2]
  #print np.fromstring(CT.PixelData, np.dtype('int16'))
  a[z] = np.fromstring(CT.PixelData, np.dtype(format_str)).reshape((N[1],N[0])) #np.dtype(format_str)
  #print a.min(), a.max()
  z+=1
a = a*int(dicom.read_file(fCT[0]).RescaleSlope)+int(dicom.read_file(fCT[0]).RescaleIntercept)
print np.min(a),np.max(a)
a = np.clip(a,-1000,2500)
b = np.zeros((N[2]+2,N[1]+2,N[0]+2),'int16')

b[1:-1,1:-1,1:-1] = a
b = b.flatten()
h.Set(len(b),b)
h.Write()
print np.min(a),np.max(a)
rhoHist = rhoFit(a)
c = np.zeros((N[2]+2,N[1]+2,N[0]+2),'float64')
c[1:-1,1:-1,1:-1] = rhoHist
c = c.flatten()
rho.Set(len(c),c)
rho.Write()
f.Close()
  
print A['DATA_PATH']+'/out/phantom.root created: check the phantom'
sys.exit()
