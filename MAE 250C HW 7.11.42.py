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


# In[3]:


# numerical integration to find cycles to failure

deltaNarray = []
b = 38.0 # mm
t = 6.0 # mm
m = 3.24
C0 = 5.11e-10 
Pmax = 120000 # N
Pmin = 40000 # N
deltaP = Pmax - Pmin
deltaS = deltaP/(2*b*t)
sig0 = 1255 # MPa
R = Pmin/Pmax
gamma = .420
C = C0/((1-R)**(m*(1-gamma)))
C = C/1000
a = np.linspace(.001,.0288,10000)
Nf = []
for i in range(0,9999):
        abar = a[i] + ((a[i+1]-a[i])/2)
        alpha = abar/.038 
        F = (1 - (0.5*alpha) + (0.326*(alpha**2)))/np.sqrt(1-alpha)
        deltaK = F*deltaS*np.sqrt(3.14159*abar)
        dadN  = C*(deltaK**m)
        deltaN = (a[i+1]-a[i])/dadN
        deltaNarray.append(deltaN)
sum(deltaNarray)


# In[ ]:




