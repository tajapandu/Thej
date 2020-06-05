#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pylab
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import Series, DataFrame
import seaborn as sns
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter
from sklearn import linear_model
import statsmodels.api as sm


# In[57]:


# numerical integration to find cycles to failure

b = 40.0 # mm
t = 20.0 # mm
m = 2.2

Smax = np.array([46.0,46.0,40.0,36.0,42.0,44.0])
Smin = np.array([0.0,20.0,24.0,24.0,18.0,28.0])
deltaS1 = Smax[0]-Smin[0]  # Takeoff
deltaS2 = Smax[1]-Smin[1] # Climb, 15 degree flaps
deltaS3 = Smax[2]-Smin[2] # Climb to cruise
deltaS4 = Smax[3]-Smin[3] # Cruise
deltaS5 = Smax[4]-Smin[4] # Descent+Hold
deltaS6 = Smax[5]-Smin[5] # Approach
R1 = 0.0/46.0
R2 = 20.0/46.0
R3 = 24.0/40.0
R4 = 24.0/36.0
R5 = 18.0/42.0
R6 = 28.0/44.0
Kc = 70.0 # ksi*sqrt(in)
C=3.0e-7
a = [.125]
c = 0.125
i = 0

def calcdeltaA(deltaS,R):
    d = a[len(a)-1]/(c + a[len(a)-1])
    F = 0.5*(3.0-d)*(1+(1.243*((1-d)**3.0)))
    deltaK = F*deltaS*np.sqrt(np.pi*a[len(a)-1])
    dadN  = (C*(deltaK**2.2))/(((1-R)*Kc)-deltaK)
    deltaA = dadN*1
    a.append(deltaA + a[len(a)-1])
    return deltaK

failure = False

while failure == False:
    
    if calcdeltaA(deltaS1,R1) >= (1-R1)*Kc:  #Takeoff (1 cycle)
        failure = True
        break
        
    for k in range(len(a)-1,len(a)+1):       #Climb, 15 degree flaps (2 cycles)
        if calcdeltaA(deltaS2,R2) >= (1-R2)*Kc:
            failure = True
            break 
    if failure == True:
        break
    
    for l in range(len(a)-1,len(a)+2):       #Climb to cruise (3 cycles)
        if calcdeltaA(deltaS3,R3) >= (1-R3)*Kc:
            failure = True
            break
    if failure == True:
        break
        
    for m in range(len(a)-1,len(a)+16):      #Cruise (17 cycles)
        if  calcdeltaA(deltaS4,R4) >= (1-R4)*Kc:
            failure = True
            break
    if failure == True:
        break
        
    for n in range(len(a)-1,len(a)+2):       #Descent + Hold (3 cycles)
        if calcdeltaA(deltaS5,R5) >= (1-R5)*Kc:
            failure = True
            break
    if failure == True:
        break
        
    for p in range(len(a)-1,len(a)+2):       #Approach (3 cycles)
        if calcdeltaA(deltaS6,R6) >= (1-R6)*Kc:
            failure = True
            break
    if failure == True:
        break
        
    i = i+1     
print(i)
print(R1,R2,R3,R4,R5,R6)


# In[ ]:





# In[19]:


(3e-7)*1000000


# In[ ]:




