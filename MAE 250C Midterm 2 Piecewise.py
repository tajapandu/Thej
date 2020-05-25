#!/usr/bin/env python
# coding: utf-8

# In[49]:


import pylab
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import Series, DataFrame
import seaborn as sns
from scipy import optimize
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter


# In[52]:


# 7075 T6 Al
# Define and input problem specifics

# Geometry
Kt = 3.00

# Material
E = 71000.0 # MPa
delSig = np.array([284.0,123.0,113.66,62.75,64.41])
delEps = np.array([0.004,0.002,0.004,0.006,0.014])
P = np.array([0.0,0.0,0.0,0.0,0.0])
sig_fp = 1466.0
b = -0.143
c = -0.619
eps_fp = 0.262
Hpr = 977.0
npr = 0.106

# Stress History

numPks = 9
numEl = 4
S = np.array([0.0,350.0,70.0,280.0,70.0,210.0,70.0,140.0,0.0,350.0])

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
            P[j] = P[j] + signS*((Sint-SL)/(delSig[j])) # update availability coefficient
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
print(sigpts)
print(epspts)


# In[59]:


deleps0 = .004
deleps1 = .006
deleps2 = .010
deleps3 = .016
deleps4 = .030
def f(delsig1):
    return ((delsig1/E) + ((delsig1/(Hpr))**(1/npr)))-deleps1
delsig1 = optimize.brentq(f,0, 1000000, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
delsig1


# In[45]:


def f(delsig2):
    return ((delsig2/E) + ((delsig2/(Hpr))**(1/npr)))-deleps2
delsig2 = optimize.brentq(f,0, 1000000, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
delsig2


# In[46]:


def f(delsig3):
    return ((delsig3/E) + ((delsig3/(Hpr))**(1/npr)))-deleps3
delsig3 = optimize.brentq(f,0, 1000000, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
delsig3


# In[47]:


def f(delsig4):
    return ((delsig4/E) + ((delsig4/(Hpr))**(1/npr)))-deleps4
delsig4 = optimize.brentq(f,0, 1000000, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
delsig4


# In[48]:


def f(delsig0):
    return ((delsig0/E) + ((delsig0/(Hpr))**(1/npr)))-deleps0
delsig0 = optimize.brentq(f,0, 1000000, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
delsig0


# In[60]:


print(sigpts)


# In[61]:


print(epspts)


# In[69]:


SWT1 = 5.17571831
SWT2 = 1.928603188
SWT3 = 0.668670287
SWT4 = 0.023771958

def f(Nf1):
    return (((sig_fp**2)/E)*((2*Nf1)**(2*b)) + (sig_fp*eps_fp*((2*Nf1)**(b+c))))-SWT1
Nf1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
Nf1


# In[68]:


def f(Nf2):
    return (((sig_fp**2)/E)*((2*Nf2)**(2*b)) + (sig_fp*eps_fp*((2*Nf2)**(b+c))))-SWT2
Nf2 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
Nf2


# In[67]:


def f(Nf3):
    return (((sig_fp**2)/E)*((2*Nf3)**(2*b)) + (sig_fp*eps_fp*((2*Nf3)**(b+c))))-SWT3
Nf3 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
Nf3


# In[66]:


def f(Nf4):
    return (((sig_fp**2)/E)*((2*Nf4)**(2*b)) + (sig_fp*eps_fp*((2*Nf4)**(b+c))))-SWT4
Nf4 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
Nf4


# In[ ]:




