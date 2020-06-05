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


# In[19]:


# numerical integration to find cycles to failure

b = 40.0 # mm
t = 20.0 # mm
m = 2.2
C0 = 3.28e-9 
sig0 = 363 # MPa
deltaS1 = 46.0 # Takeoff
deltaS2 = 26.0 # Climb, 15 degree flaps
deltaS3 = 16.0 # Climb to cruise
deltaS4 = 12.0 # Cruise
deltaS5 = 24.0 # Descent+Hold
deltaS6 = 16.0 # Approach
R1 = 0.0
R2 = 0.43
R3 = 0.60
R4 = 0.67
R5 = 0.43
R6 = 0.64
Kc = 70 # ksi*sqrt(in)
C=3.0e-7
a = [.125]
c = 0.125
i = 0

def calcdeltaA(deltaS,R):
    d = a[len(a)-1]/(c + a[len(a)-1])
    F = 0.5*(3.0-d)*(1+(1.243*((1-d)**3)))
    deltaK = F*deltaS*np.sqrt(3.14159*a[len(a)-1])
    dadN1  = (C*(deltaK**2.2))/(((1-R)*Kc)-deltaK)
    deltaA1 = dadN1*1
    a.append(deltaA1 + a[len(a)-1])
    return deltaK

failure = False

while failure == False:
    if calcdeltaA(deltaS1,R1) >= (1-R1)*Kc:
        failure = True
        break
    
    for k in range(len(a)-1,len(a)+1): 
        if calcdeltaA(deltaS2,R2) >= (1-R2)*Kc:
            failure = True
            break 
    if failure == True:
        break
    
    for l in range(len(a)-1,len(a)+2):
        if calcdeltaA(deltaS3,R3) >= (1-R3)*Kc:
            failure = True
            break
    if failure == True:
        break
        
    for m in range(len(a)-1,len(a)+16):
        if  calcdeltaA(deltaS4,R4) >= (1-R4)*Kc:
            failure = True
            break
    if failure == True:
        break
        
    for n in range(len(a)-1,len(a)+2):
        if calcdeltaA(deltaS5,R5) >= (1-R5)*Kc:
            failure = True
            break
    if failure == True:
        break
        
    for p in range(len(a)-1,len(a)+2):
        if calcdeltaA(deltaS6,R6) >= (1-R6)*Kc:
            failure = True
            break
    if failure == True:
        break
        
    i = i+1     
print(i)


# In[ ]:





# In[19]:


(3e-7)*1000000


# In[ ]:




