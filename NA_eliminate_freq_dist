#!/usr/bin/jupyterenv python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cooler
import matplotlib.colors


#open and read mcool file with selectable resolution
c = cooler.Cooler('1mb_mm9_woSc_B.mcool'"::/resolutions/1000000") ## selection of resolution

pixs = c.matrix(as_pixels=True).fetch('chr1') ##selecting the chr -- make it user input at beginning
GAM = pd.DataFrame(data=pixs)
pix_weight = c.bins().fetch('chr1')    ##selecting the chr -- make it user input at beginning
temp = pix_weight.loc[:,'weight']
qGAM = len (temp) - 1

##create a empty array
tempGAM = pixs
pGAM = len (tempGAM)
tGAM = tempGAM.iloc[pGAM-1,1] - tempGAM.iloc[0,0]
tGAM = int(tGAM)
index=range(0,qGAM)
columns=range(0,6)
emptyGAM = pd.DataFrame(index=index, columns=columns)

for i in range (0,qGAM):
    emptyGAM.iloc[i,0] = i+1      #difference
    emptyGAM.iloc[i,1] = 0        #freq(sum) of difference unbalance
    emptyGAM.iloc[i,2] = qGAM-i     #number of total element
    emptyGAM.iloc[i,3] = 0        #freq mean  unbalance
    emptyGAM.iloc[i,4] = 0        #freq(sum) of difference balance
    emptyGAM.iloc[i,5] = 0        #freq mean  balance

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
    
##eliminate NA values 
for i in range (0,qGAM):
    if pd.isnull(pix_weight.loc[i,'weight'])== True:
        emptyGAM.iloc[0:qGAM-i,2] -= 1
        
emptyGAM = emptyGAM[emptyGAM.iloc[:,2] >= 1]  ##here is the selection for min difference -- make it initially selected

##calculate mean values -- if not values will be too high -- not must but suggested
emptyGAM.iloc[:,3] = emptyGAM.iloc[:,1]/emptyGAM.iloc[:,2]
emptyGAM.iloc[:,5] = emptyGAM.iloc[:,4]/emptyGAM.iloc[:,2]

"""
#writing the result on a csv file 
meanGAM =open("results/RESULTS_freq_1mb_chr1_1mRes_0313.csv" , "a")

for i in range (0, tGAM):
    meanGAM.write (str(emptyGAM.iloc[i,1])) ##unbalanced sum
    meanGAM.write ('\t')
    meanGAM.write (str(emptyGAM.iloc[i,3])) ##unbalanced mean
    meanGAM.write ('\t')
    meanGAM.write (str(emptyGAM.iloc[i,4])) ##balanced sum
    meanGAM.write ('\t')
    meanGAM.write (str(emptyGAM.iloc[i,5])) ##balanced mean
    meanGAM.write ('\n')
    
meanGAM.close()
 """   
##same code with upper part 
##change the file save location with aotomated one -- like give and auto generated name based on file name, chr number , date etc
RESULT = pd.DataFrame(data=emptyGAM.iloc[:,[1,3,4,5]])
RESULT.to_csv('results/RESULTS_freq_1mb_chr1_1mRes_0313-1.csv', index=False , header=None , sep='\t')    


#plotting the freq-distance graph
x = emptyGAM.iloc[1:-1,3]
y = emptyGAM.iloc[1:-1,5]*50000  ## to many difference between balance-unbalance data so need to multiply make user enter

plt.switch_backend('agg')
plt.figure(figsize=(20,5))
plt.plot(x.iloc[20:,]  ,color='green', label='unbalanced')
plt.plot(y.iloc[20:,] ,color='blue', label='balanced')
plt.legend(loc='upper right')

plt.xlabel('distance difference')
plt.ylabel('frequency')
plt.savefig('results/NA_eliminated_freq_1mb_chr1_1mRes_0313.png')
