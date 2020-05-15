#!/usr/bin/env python
# coding: utf-8

# In[76]:


import pylab
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import Series, DataFrame
import seaborn as sns
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter


# In[ ]:





# In[263]:


# 7075 T6 Al
# Define and input problem specifics

# Geometry
Kt = 4.00

# Material
E = 71000.0 # MPa
delSig = np.array([400,125,50,50,75])
delEps = np.array([0.00563,0.00437,0.00480,0.00920,0.02900])
P = np.array([0.0,0.0,0.0,0.0,0.0])

# Stress History

numPks = 19
numEl = 4
S = np.array([0,336,48.3,292,48.3,255,48.3,219,48.3,181,48.3,143,48.3,105,48.3,67.9,-19.3,48.3,-67.9,336])

# Set Initial Conditions

E0 = S0 = 0
EL = SL = 0
signS = 1
sigpts = []
epspts = []

# Loop over turning points

for i in range (1,numPks+1):
    E0 = E0 + signS*EL                           #Set origin for this stress range
    S0 = S0 + signS*SL
    EL = 0                                       #Set initial local coordinates to (0,0) for this turning point                               
    SL = 0
    delS = (S[i] - S[i-1])                       #Change of stress between turning points
    signS = delS/(abs(delS))                     #Sign of stress between turning points (either +1 or -1)
    sigEps = ((Kt*delS)**2)/E                    #Neuber Hyperbola
    for j in range (0,numEl+1):
        Str = SL + abs(signS - P[j])*delSig[j]   #Stretch element to max
        Etr = EL + abs(signS - P[j])*delEps[j]
        sigEpsTr = Str*Etr
        if sigEpsTr < sigEps:
            P[j] = signS
            EtrGlobal = E0 + signS*Etr
            StrGlobal = S0 + signS*Str
            sigpts.append(StrGlobal)
            epspts.append(EtrGlobal)
            SL = Str
            EL = Etr
        else:
            m = (Str-SL)/(Etr-EL)
            Eint = (((m*EL)-SL)+(np.sqrt(((SL-(m*EL))**2)+(4*m*sigEps))))/(2*m) # Find hyperbola/element intersection
            Sint = SL + m*(Eint-EL)
            P[j] = P[j] + signS*((Sint-SL)/(delSig[j]))                         # update availability coefficient
            SintGlobal = S0 + signS*Sint                                        # store in global coordinates
            EintGlobal = E0 + signS*Eint
            epspts.append(EintGlobal)                                           # add global coordinates to an array
            sigpts.append(SintGlobal)
            SL = Sint
            EL = Eint
            break

sigpts = np.insert(sigpts,0,0)
epspts = np.insert(epspts,0,0)
sns.set_style('whitegrid')
plt.plot(epspts,sigpts,'g') 
plt.ylabel(r'$\sigma$')
plt.xlabel(r'$\epsilon$')
plt.show()

