# -*- coding: utf-8 -*-
"""
Created on Fri Sep 29 18:04:17 2023

@author: mikej
"""

"""
Not going to bed tongiht until:
    1. finish part (a)
    2. have a good idea and have started (b)
    3. have an understanding of how to do the rest of the parts
"""

import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
from scipy.constants import c

def uvscatterplot(xdata,ydata,plot):
    if plot == True:
        plt.figure(figsize = (8,8))
        plt.scatter(xdata, ydata, s = 50, marker = '.')
        plt.axhline(0, color = 'gray', linestyle = '--')
        plt.axvline(0, color = 'gray', linestyle = '--')
        plt.xlabel("u ($\lambda$)", fontsize = 20)
        plt.ylabel("v ($\lambda$)", fontsize = 20)
        plt.grid(True)
    else:
        pass

def DataKeyMapper(key):
    if key == 'a':
        Adata = np.loadtxt("vla_a_array.txt")
    elif key == 'd':
        Adata = np.loadtxt("vla_d_array.txt")
    return Adata

Adata = DataKeyMapper('a')
#a_data1 = pd.read_csv("vla_a_array.txt", delimiter= '\t') #if want to load in using pandas
   
Adata = Adata[:,:3]  #only the first three columns contain coords
Adata = Adata*1e-9*c
Nant = Adata.shape[0]
Nvis = Nant*(Nant - 1)//2 #number of visibilities = NC2 -> floor function incase Nant is odd
print(Nvis)
Lat = 34.00784*np.pi/180

zenith = np.array([np.cos(Lat), 0, np.sin(Lat)])
east = np.array([0,1,0])
north = np.cross(zenith, east)


"""
Will need to particularly write this section below out by hand to make 
sure I truly understand the matrix multiplication
"""

#coords in the following definition is of unit length in each direction
coords = np.vstack([north, east, zenith])   #define the coordinate system in NS/EW coords

xyz = Adata[:,:3]@coords
if True:
    plt.figure(figsize = (8,8))
    plt.scatter(xyz[:,0],xyz[:,1], s = 125, marker = 'x')
    plt.axhline(0, color = 'gray', linestyle = '--')
    plt.axvline(0, color = 'gray', linestyle = '--')
    plt.grid(True)
    


"Print out the RMS error or vertical scatter"
print("The spread in the E/W direction is", xyz[:,0].std())
print("Spread in the N/S direction is", xyz[:,1].std())

z = xyz[:,2]
Vscat = np.sum((z - z.mean())**2)/len(z)

print("vertical RMS scatter is", Vscat)



"plot out the UV plots"
uv = np.zeros([Nvis,2])

start = time.time()
istep = 0
for i in range(Nant):
    for j in range(i+1,Nant):
        uv[istep,:] = xyz[i,:2] - xyz[j,:2]
        istep += 1
uv = np.vstack([uv,-1*uv])
stop = time.time()
# print(stop - start)  #benchmarking vectorized version
print(uv.shape)

lmda = c/1.4e9
uv[:,0], uv[:,1] = uv[:,0]*100/21.4, uv[:,1]*100/21.4
uvscatterplot(uv[:,0], uv[:,1], True)




































