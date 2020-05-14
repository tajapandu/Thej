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





# In[181]:


# 7075 T6 Al
# Define and input problem specifics

# Geometry
Kt = 4.00
w1 = .03810 # m
w2 = .05715 # m
rho = .00145 # m
t = .00229 # m

# Material
E = 71000.0 # MPa
delSig = np.array([400,125,50,50,75])
delEps = np.array([0.00563,0.00437,0.00480,0.00920,0.02900])
P = np.array([0,0,0,0,0])

# Stress History

numPks = 19
numEl = 4
S = np.array([0,336,48.3,292,48.3,255,48.3,219,48.3,181,48.3,143,48.3,105,48.3,67.9,-19.3,48.3,-67.9,336])

# Set Initial Conditions

E0 = S0 = 0
EL = SL = 0
signS = 1

# Loop over turning points

for i in range (1,numPks+1):
    E0 = E0 + signS*EL          #Set origin for this stress range
    S0 = S0 + signS*SL
    EL = 0                      #Set initial local coordinates to (0,0) for this turning point                               
    SL = 0
    delS = (S[i] - S[i-1]) #Change of stress between turning points
    signS = delS/(abs(delS)) #Sign of stress between turning points (either +1 or -1)
    sigEps = ((Kt*delS)**2)/E   #Neuber Hyperbola
    for j in range (0,numEl+1):
        Str = SL + abs(signS - P[j])*delSig[j]  #Stretch element to max
        Etr = EL + abs(signS - P[j])*delEps[j]
        sigEpsTr = Str*Etr
        if sigEpsTr < sigEps:
            P[j] = signS
            EtrGlobal = E0 + signS*Etr
            StrGlobal = S0 + signS*Str
            #plot segment from (EL,SL) to (Etr,Str) IN GLOBAL COORDINATES
            #plt.plot ([EL,EtrGlobal],[SL,StrGlobal])
            SL = Str
            EL = Etr
        else:
            m = (Str-SL)/(Etr-EL)
            Eint = (((m*EL)-SL)+(np.sqrt(((SL-(m*EL))**2)+(4*m*sigEps))))/(2*m)
            Sint = SL + m*(Eint-EL)
            P[j] = P[j] + signS*((Sint-SL)/(delSig[j]))
            #write or store (Sint,Eint) IN GLOBAL COORDINATES
            SintGlobal = S0 + signS*Sint
            EintGlobal = E0 + signS*Eint
            #plot segment from (EL,SL) to (Eint,Sint) in global coordinates
            #lt.plot([EL,EintGlobal],[SL,SintGlobal])
            SL = Sint
            EL = Eint
            break
            
            


# In[134]:





# In[113]:





# In[150]:


df = pd.read_excel(r'C:\Users\tajap\Desktop\HW4_GES.xlsx')


# In[151]:


df


# In[176]:


plt.plot([df.E[0],df.E[1]],[df.S[0],df.S[1]])
plt.plot([df.E[1],df.E[2]],[df.S[1],df.S[2]])
plt.plot([df.E[2],df.E[3]],[df.S[2],df.S[3]])
plt.plot([df.E[3],df.E[4]],[df.S[3],df.S[4]])
plt.plot([df.E[4],df.E[5]],[df.S[4],df.S[5]])
plt.plot([df.E[5],df.E[6]],[df.S[5],df.S[6]])
plt.plot([df.E[6],df.E[7]],[df.S[6],df.S[7]])
plt.plot([df.E[7],df.E[8]],[df.S[7],df.S[8]])
plt.plot([df.E[8],df.E[9]],[df.S[8],df.S[9]])
plt.plot([df.E[9],df.E[10]],[df.S[9],df.S[10]])
plt.plot([df.E[10],df.E[11]],[df.S[10],df.S[11]])
plt.plot([df.E[11],df.E[12]],[df.S[11],df.S[12]])
plt.plot([df.E[12],df.E[13]],[df.S[12],df.S[13]])
plt.plot([df.E[13],df.E[14]],[df.S[13],df.S[14]])
plt.plot([df.E[14],df.E[15]],[df.S[14],df.S[15]])
plt.plot([df.E[15],df.E[16]],[df.S[15],df.S[16]])
plt.plot([df.E[16],df.E[17]],[df.S[16],df.S[17]])
plt.plot([df.E[17],df.E[18]],[df.S[17],df.S[18]])
plt.plot([df.E[18],df.E[19]],[df.S[18],df.S[19]])
plt.plot([df.E[19],df.E[20]],[df.S[19],df.S[20]])
plt.plot([df.E[20],df.E[21]],[df.S[20],df.S[21]])
plt.plot([df.E[21],df.E[22]],[df.S[21],df.S[22]])
plt.plot([df.E[22],df.E[23]],[df.S[22],df.S[23]])
plt.plot([df.E[23],df.E[24]],[df.S[23],df.S[24]])
plt.plot([df.E[24],df.E[25]],[df.S[24],df.S[25]])
plt.plot([df.E[25],df.E[26]],[df.S[25],df.S[26]])
plt.plot([df.E[26],df.E[27]],[df.S[26],df.S[27]])
plt.plot([df.E[27],df.E[28]],[df.S[27],df.S[28]])
plt.plot([df.E[28],df.E[29]],[df.S[28],df.S[29]])
plt.plot([df.E[29],df.E[30]],[df.S[29],df.S[30]])
plt.plot([df.E[30],df.E[31]],[df.S[30],df.S[31]])
plt.plot([df.E[31],df.E[32]],[df.S[31],df.S[32]])
plt.plot([df.E[32],df.E[33]],[df.S[32],df.S[33]])
plt.plot([df.E[33],df.E[34]],[df.S[33],df.S[34]])
plt.plot([df.E[34],df.E[35]],[df.S[34],df.S[35]])
plt.plot([df.E[35],df.E[36]],[df.S[35],df.S[36]])
plt.plot([df.E[36],df.E[37]],[df.S[36],df.S[37]])
plt.plot([df.E[37],df.E[38]],[df.S[37],df.S[38]])


# In[ ]:


for n in df:
    plt.plot([df.Global_Strain[n],df.GlobalStress[n]],[df.Global_Strain[n+1],df.Global_Stress[n+1]])

