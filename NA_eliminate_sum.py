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


#read mcool file
c = cooler.Cooler('/../../n/scratch2/onur/MARCH/4DNFI7Z381GX.mcool'"::/resolutions/100000") #need to user input bot file and resolution


pixs = c.matrix(as_pixels=True).fetch('chr19')
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
meanGAM =open("RESULTS_sum_chr19_4DN_0308.txt" , "a")

for i in range (1, tGAM-1):
    meanGAM.write (str(emptyGAM.iloc[i,1])) ##unbalanced sum
    meanGAM.write ('\t')
    meanGAM.write (str(emptyGAM.iloc[i,3])) ##unbalanced mean
    meanGAM.write ('\t')
    meanGAM.write (str(emptyGAM.iloc[i,4])) ##balanced sum
    meanGAM.write ('\t')
    meanGAM.write (str(emptyGAM.iloc[i,5])) ##balanced mean
    meanGAM.write ('\n')
    
meanGAM.close()


##plotting

x = emptyGAM.iloc[:,1]
y = emptyGAM.iloc[:,4]*10000   ##balance is too small to observe

#plt.yscale('log')
plt.figure(figsize=(15,4))
plt.plot(x.iloc[20:-1]  ,color='green', label='unbalanced') 
plt.plot(y.iloc[20:-1] ,color='blue', label='balanced')
plt.legend(loc='upper right')

plt.xlabel('location')
plt.ylabel('sum')
#plt.show()
plt.savefig('NA-eliminated_sum_chr19_4DN_0308.png')
