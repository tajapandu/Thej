#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[15]:


# numerical integration to find cycles to failure

K1array = []
K2array = []
K3array = []
deltaAarray = []
b = 40.0 # mm
t = 20.0 # mm
m = 3.13
C0 = 3.28e-9 
sig0 = 363 # MPa
ai = .001
af = .0143
deltaS1 = 112.5 # MPa
deltaS2 = 56.25 # MPa
deltaS3 = 225 # MPa
R1 = 0.5
R2 = 0
R3 = 0
gamma = .928
C1 = C0/((1-R1)**(m*(1-gamma)))
C1 = C1/1000
C2 = C0/((1-R2)**(m*(1-gamma)))
C2 = C2/1000
C3 = C0/((1-R3)**(m*(1-gamma)))
C3 = C3/1000
F = 1.12
a = [.001]
i = 0
while a[i] <= 0.014:
    deltaK1 = F*deltaS1*np.sqrt(3.14159*a[i])
    dadN1  = C1*(deltaK1**m)
    deltaA1 = dadN1*200
      
   
    deltaK2 = F*deltaS2*np.sqrt(3.14159*a[i])
    dadN2  = C2*(deltaK2**m)
    deltaA2 = dadN2*1000
    
    deltaK3 = F*deltaS3*np.sqrt(3.14159*a[i])
    dadN3  = C3*(deltaK3**m)
    deltaA3 = dadN3*1
    
    deltaA = deltaA1 + deltaA2 + deltaA3
    
    a.append(deltaA+a[i])
        
    i = i+1
print(i)


# In[28]:


print(C1,C2,C3)


# In[ ]:


dadB = sum(deltaAarray)
Blocks = (af-ai)/dadB
Blocks

