#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cooler
import matplotlib.colors
from itertools import combinations
import fileinput
import os
import sys
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

#open and read mcool file with selectable resolution
c = cooler.Cooler('%s'"::/resolutions/%s" %(file_name,res)) 
pixs = c.matrix(as_pixels=True).fetch('chr%s' %limit)
pix_weight = c.bins().fetch('chr%s' %limit)
temp = pix_weight.loc[:,'weight']
qGAM = len (temp) - 1

##create a empty array
tempGAM = pixs
pGAM = len (tempGAM)

##create a empty array
tGAM = tempGAM.iloc[pGAM-1,1] - tempGAM.iloc[0,0]
tGAM = int(tGAM)
index=range(0,qGAM)
columns=range(0,6)
emptyGAM = pd.DataFrame(index=index, columns=columns)
for i in range (0,qGAM):
    emptyGAM.iloc[i,0] = i+1      #difference
    emptyGAM.iloc[i,1] = 0        #freq(sum) of difference unbalance
    emptyGAM.iloc[i,2] = qGAM-i   #number of total element
    emptyGAM.iloc[i,3] = 0        #freq mean  unbalance
    emptyGAM.iloc[i,4] = 0        #freq(sum) of difference balance
    emptyGAM.iloc[i,5] = 0        #freq mean

##sum and mean of the distance values include 0

for i in range (0,pGAM):
    lGAM = tempGAM.iloc[i,1] - tempGAM.iloc[i,0] -1
    lGAM = int(lGAM)
    
    if pd.isnull(tempGAM.iloc[i,3])== False: 
        tB = tempGAM.iloc[i,3]
        emptyGAM.iloc[lGAM,1] += tempGAM.iloc[i,2]
        
    if pd.isnull(tempGAM.iloc[i,3])== True: 
        tB = 0.0
          
    tB = float(tB)
    emptyGAM.iloc[lGAM,4] += tB

    
#find weigth column
def column_index(df, query_cols):
    cols = df.columns.values
    sidx = np.argsort(cols)
    return sidx[np.searchsorted(cols,query_cols,sorter=sidx)]

a = column_index(pix_weight, ['weight'])
a = int(a)
    
##eliminate NA values 
for i in range (0,qGAM):
    if pd.isnull(pix_weight.loc[i,'weight'])== True:
        emptyGAM.iloc[0:qGAM-i,2] -= 1
        
emptyGAM = emptyGAM[emptyGAM.iloc[:,2] >= 1]  ##here is the selection for min difference -- make it initially selected

##calculate mean values -- if not values will be too high -- not must but suggested
emptyGAM.iloc[:,3] = emptyGAM.iloc[:,1]/emptyGAM.iloc[:,2]
emptyGAM.iloc[:,5] = emptyGAM.iloc[:,4]/emptyGAM.iloc[:,2]


##change the file save location with aotomated one -- like give and auto generated name based on file name, chr number , date etc
RESULT = pd.DataFrame(data=emptyGAM.iloc[:,[1,3,4,5]])
RESULT.to_csv('results/RESULTS_freq_1mb_chr1_1mRes_0313-1.csv', index=False , header=None , sep='\t')    


#plotting the freq-distance graph
x = emptyGAM.iloc[1:-1,3]
y = emptyGAM.iloc[1:-1,5]*500 ##this can be changed

plt.switch_backend('agg')
plt.figure(figsize=(20,5))
plt.plot(x.iloc[1:,]  ,color='green', label='unbalanced')
plt.plot(y.iloc[1:,] ,color='blue', label='balanced')
plt.legend(loc='upper right')

plt.xlabel('distance difference')
plt.ylabel('frequency')
plt.savefig('%s_BxuB.png' %out_name)

##for seperate plots 
x = emptyGAM.iloc[1:-1,3]
y = emptyGAM.iloc[1:-1,5]

#unbalanced
plt.switch_backend('agg')
plt.figure(figsize=(20,5))
plt.plot(x.iloc[1:,]  ,color='green', label='unB_chr%s_800k_COV' %limit)
plt.legend(loc='upper right')
plt.xlabel('distance difference')
plt.ylabel('frequency')
plt.savefig('%s_unBal.png' %out_name)

#balanced
plt.switch_backend('agg')
plt.figure(figsize=(20,5))
plt.plot(y.iloc[1:,] ,color='blue', label='B_chr%s_800k_COV' %limit)
plt.legend(loc='upper right')
plt.xlabel('distance difference')
plt.ylabel('frequency')
plt.savefig('%s_Bal.png' %out_name)
