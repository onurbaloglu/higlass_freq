#!/usr/bin/env python

import pandas as pd
import numpy as np
import math
import statistics
import matplotlib.pyplot as plt
import os
import sys
import cooler
import matplotlib.colors
from itertools import combinations
import fileinput
import os.path

for i in range (1,5):
    file_name = input ("enter the file location and name of the mcool file:\n")
    if os.path.isfile(file_name):
        print ('your file name is: %s' %file_name)
        break
    else :
        print ("folder not found!! enter the folder name of the mcool file again:\n  ")
        if i == 4:
            print('iteration exceeded!!! ')
limit = int(input('enter the chr number:\n'))
out_name =str(input("enter the basename for the output file:\n"))
res = int(input('enter the resolution :\n'))

#read mcool file
c = cooler.Cooler('%s'"::/resolutions/%s" %(file_name,res)) 

pixs = c.matrix(as_pixels=True).fetch('chr%s' %limit)
sub = pixs.iloc[0,0]
pixs.iloc[:,0] = pixs.iloc[:,0] - sub
pixs.iloc[:,1] = pixs.iloc[:,1] - sub
pix = pixs.iloc[:]



###opening empty file

tempGAM = pixs
qGAM = len (tempGAM.columns) #column size
pGAM = len (tempGAM)        #raw size
tGAM = tempGAM.iloc[pGAM-1,1] - tempGAM.iloc[0,0]
tGAM = int(tGAM)
index=range(0,tGAM+1)
columns=range(0,6)
emptyGAM = pd.DataFrame(index=index, columns=columns)

for i in range (0,tGAM+1):
    emptyGAM.iloc[i,0] = i+1  #bin id
    emptyGAM.iloc[i,1] = 0    #sum of unbalanced
    emptyGAM.iloc[i,2] = tGAM-1 #number of elements-  '-1' means we start from difference 2 make it user input
    emptyGAM.iloc[i,3] = 0    #mean unblance
    emptyGAM.iloc[i,4] = 0    #sum of balanced
    emptyGAM.iloc[i,5] = 0    #mean blance
    
##writing sums on file    
for i in range (0,pGAM):
    x = int(tempGAM.iloc[i,0])
    y = int(tempGAM.iloc[i,1])
    if y-x >= 2:
        tuB = tempGAM.iloc[i,2]
        tuB = float(tuB)
        
        emptyGAM.iloc[x,1] += tuB
        emptyGAM.iloc[y,1] += tuB
        
        if pd.isnull(tempGAM.iloc[i,3])== False:
            tB = tempGAM.iloc[i,3]
    
        if pd.isnull(tempGAM.iloc[i,3])== True:
            tB = 0.0
        tB = float(tB)
        emptyGAM.iloc[x,4] += tB
        emptyGAM.iloc[y,4] += tB

emptyGAM.iloc[:,3] = emptyGAM.iloc[:,1]/emptyGAM.iloc[:,2]
emptyGAM.iloc[:,5] = emptyGAM.iloc[:,4]/emptyGAM.iloc[:,2]


##writing a result file
RESULT = pd.DataFrame(data=emptyGAM.iloc[:,[1,3,4,5]])
RESULT.to_csv('%s.csv' %out_name, index=False , header=None , sep='\t' )


##plotting
##for seperate plots 
x = emptyGAM.iloc[1:-1,3]
y = emptyGAM.iloc[1:-1,5]

#unbalanced
plt.switch_backend('agg')
plt.figure(figsize=(20,5))
plt.plot(x.iloc[5:,]  ,color='green', label= '%s_unBal.png'  %out_name)
plt.legend(loc='upper right')
plt.xlabel('location')
plt.ylabel('sum')
plt.savefig('%s_unBal.png'  %out_name)

#balanced
plt.switch_backend('agg')
plt.figure(figsize=(20,5))
plt.plot(y.iloc[5:,] ,color='blue', label= '%s_Bal.png'  %out_name)
plt.legend(loc='upper right')
plt.xlabel('location')
plt.ylabel('sum')
plt.savefig('%s_Bal.png'  %out_name)

##plotting together

x = emptyGAM.iloc[1:-1,3]
y = emptyGAM.iloc[1:-1,5]*1000

#plt.yscale('log')
plt.figure(figsize=(15,4))
plt.title('%s' %out_name, size=18)
plt.plot(x.iloc[20:-1]  ,color='green', label='unbalanced')
plt.plot(y.iloc[20:-1] ,color='blue', label='balanced')
plt.legend(loc='upper right')
plt.xticks(size=12)
plt.yticks(size=12)

plt.xlabel('location')
plt.ylabel('sum')
#plt.show()
plt.savefig('%s_BuB.png'  %out_name)
