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
from scipy import optimize
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter


# In[549]:


# Geometry

Kt = 4.00

# Material

E = 71000.0
sig_fp = 1466.0
b = -0.143
c = -0.619
eps_fp = 0.262
Hpr = 977.0
npr = 0.106


# Stress History

numPks = 19
S = np.array([0.0,336.0,48.3,292.0,48.3,255.0,48.3,219.0,48.3,181.0,48.3,143.0,48.3,105.0,48.3,67.9,-19.3,48.3,-67.9,336.0])
SigOr = np.array([0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0])

# Data Storage Arrays
sigPk = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
epsPk = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
delsig = [0.0]

# Handle first cycle

delS = S[1]-S[0]
signS = delS/(abs(delS))
delS = abs(delS)
# find delsig using non-linear root finding
def f(sigma):
        return ((1/Kt)*np.sqrt((sigma**2)+(sigma*E*((sigma/Hpr)**(1/npr))))-delS)
delsig = optimize.brentq(f,delS, Kt*delS, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
delsig = delsig[0]
deleps = (delsig/E)+((delsig/Hpr)**(1/npr))
sigPk[1] = sigPk[0] + signS*delsig
epsPk[1] = epsPk[0] + signS*deleps

# Loop over subsequent turning points 

for i in range (2,numPks+1):
    if SigOr[i] == 0:
        iOrigin = i-1
    else:
        iOrigin = SigOr[i]
    delS = S[i] - S[iOrigin]
    signS = delS/(abs(delS))
    delS = abs(delS)
    def f(sigma1):
        return ((1/Kt)*np.sqrt(((sigma1/2)**2)+((sigma1*E)/2)*((sigma1/(2*Hpr))**(1/npr)))-(delS/2))
    delsig = optimize.brentq(f,delS, Kt*delS, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
    delsig = delsig[0]
    deleps = 2*((delsig/(2*E))+((delsig/(2*Hpr))**(1/npr)))
    sigPk[i] = sigPk[iOrigin] + signS*delsig
    epsPk[i] = epsPk[iOrigin] + signS*deleps
    
# Set up a plot

Npts = 50

# plot initial loading

stress = np.linspace(sigPk[0],sigPk[1],num = Npts)
strain = (stress/E) + ((stress/Hpr)**(1/npr))
plt.plot(strain,stress,'r',linestyle='dashed')
plt.scatter(strain[49],stress[49],marker = 'x')
    
# loop over subsequent turning points

for i in range(2,numPks+1):
    if SigOr[i] == 0:
        iOrigin = i - 1
    else:
        iOrigin = SigOr[i]
    delsig = sigPk[i] - sigPk[iOrigin]
    signS = delsig/(abs(delsig))
    delsig = abs(delsig)
    stress = np.linspace(0,delsig,Npts)
    strain = (stress/E) + 2*((stress/2/Hpr)**(1/npr))
    strain = epsPk[iOrigin] + (signS*strain)
    stress = sigPk[iOrigin] + (signS*stress)
    for j in range(0,Npts):
        plt.plot(strain,stress)


sns.set_style('whitegrid')
plt.xlabel(r'$\epsilon$')
plt.ylabel(r'$\sigma$')
plt.scatter(epsPk,sigPk,marker = 'x')
plt.show()


# In[437]:


# Life Calculation using Modified Morrow equation
eps_array =[]
Damage = []
for k in range(0,numPks):
    sig_m = (sigPk[k]+sigPk[k+1])/2
    #print(sig_m)
    eps_a = (epsPk[k+1]-(epsPk[k]))/2
    eps_a = abs(eps_a)
    def f(Nf):
        return (((sig_fp/E)*(1-(sig_m/sig_fp))*((2*Nf)**b))+(eps_fp*((2*Nf)**c)))-eps_a
    Nf = optimize.brentq(f,1, 1e12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
    Nf = Nf[0]
    D = (1/Nf)
    Damage.append(D)
TotalDamage = sum(Damage)
Life = 1/TotalDamage

print(Life)

