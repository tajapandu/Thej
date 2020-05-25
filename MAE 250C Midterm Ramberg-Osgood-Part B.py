#!/usr/bin/env python
# coding: utf-8

# In[10]:


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


# In[11]:


# Geometry

Kt = 3.00

# Material 

E = 71000.0
sig_fp = 1466.0
b = -0.143
c = -0.619
eps_fp = 0.262
Hpr = 977.0
npr = 0.106


# Stress History

numPks = 9
S = np.array([0.0,350.0,70.0,280.0,70.0,210.0,70.0,140.0,0.0,350.0])
SigOr = np.array([0,0,0,0,1,0,1,0,1,0])

# Data Storage Arrays
sigPk = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
epsPk = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
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


# In[12]:


# Tensile Proof Loads
ProofS = [350,385,420,455,483,525]


# In[46]:


# Establish new Stress/Strain Origin after proof loading 100%

# Nonlinear root find for sigma of first load
def f(PLsig1):
    return ((1/Kt)*np.sqrt((PLsig1**2)+(PLsig1*E*((PLsig1/Hpr)**(1/npr))))-ProofS[0])
PLdelsig1 = optimize.brentq(f,ProofS[0], Kt*ProofS[0], args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
PLdelsig1

# Nonlinear root find for sigma of unload (factor of two expansion)
def f(PUsig1):
    return ((1/Kt)*np.sqrt(((PUsig1/2)**2)+((PUsig1*E)/2)*((PUsig1/(2*Hpr))**(1/npr)))-(ProofS[0]/2))
PUdelsig1 = optimize.brentq(f,ProofS[0], Kt*ProofS[0], args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
PUdelsig1

# Find sigma at new origin by subtracting unload stress from load stress
ProofSigma1 = PLdelsig1[0] - PUdelsig1[0]

# Find Epsilon at new origin with R.O. equations for loading and unloading
PLeps = (PLdelsig1[0]/E)+((PLdelsig1[0]/Hpr)**(1/npr))
PUeps = (PUdelsig1[0]/E) + 2*((PUdelsig1[0]/2/Hpr)**(1/npr))

# Subtract unload strain from load strain to find strain at new origin
ProofEpsilon1 = PLeps - PUeps

# New sigmax*eps_a values based on adjusting origin of S/S curve
SWTP1 = 2.378514474
SWTP2 = 0.458852466

# Nonlinear root find for Nf for each cycle
def f(NfP1):
    return (((sig_fp**2)/E)*((2*NfP1)**(2*b)) + (sig_fp*eps_fp*((2*NfP1)**(b+c))))-SWTP1
NfP1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)


def f(NfP2):
    return (((sig_fp**2)/E)*((2*NfP2)**(2*b)) + (sig_fp*eps_fp*((2*NfP2)**(b+c))))-SWTP2
NfP2 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

print(ProofEpsilon1)
print(ProofSigma1)
print(NfP1[0],NfP2[0])


# In[47]:


# Establish new Stress/Strain Origin after proof loading 110%

# Nonlinear root find for sigma of first load
def f(PLsig1):
    return ((1/Kt)*np.sqrt((PLsig1**2)+(PLsig1*E*((PLsig1/Hpr)**(1/npr))))-ProofS[1])
PLdelsig1 = optimize.brentq(f,ProofS[1], Kt*ProofS[1], args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
PLdelsig1

# Nonlinear root find for sigma of unload (factor of two expansion)
def f(PUsig1):
    return ((1/Kt)*np.sqrt(((PUsig1/2)**2)+((PUsig1*E)/2)*((PUsig1/(2*Hpr))**(1/npr)))-(ProofS[1]/2))
PUdelsig1 = optimize.brentq(f,ProofS[1], Kt*ProofS[1], args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
PUdelsig1

# Find sigma at new origin by subtracting unload stress from load stress
ProofSigma1 = PLdelsig1[0] - PUdelsig1[0]

# Find Epsilon at new origin with R.O. equations for loading and unloading
PLeps = (PLdelsig1[0]/E)+((PLdelsig1[0]/Hpr)**(1/npr))
PUeps = (PUdelsig1[0]/E) + 2*((PUdelsig1[0]/2/Hpr)**(1/npr))

# Subtract unload strain from load strain to find strain at new origin
ProofEpsilon1 = PLeps - PUeps

# New sigmax*eps_a values based on adjusting origin of S/S curve
SWTP1 = 2.085719623
SWTP2 = 0.297208142

# Nonlinear root find for Nf for each cycle
def f(NfP1):
    return (((sig_fp**2)/E)*((2*NfP1)**(2*b)) + (sig_fp*eps_fp*((2*NfP1)**(b+c))))-SWTP1
NfP1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)


def f(NfP2):
    return (((sig_fp**2)/E)*((2*NfP2)**(2*b)) + (sig_fp*eps_fp*((2*NfP2)**(b+c))))-SWTP2
NfP2 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

print(ProofEpsilon1)
print(ProofSigma1)
print(NfP1[0],NfP2[0])


# In[51]:


# Establish new Stress/Strain Origin after proof loading 120%

# Nonlinear root find for sigma of first load
def f(PLsig1):
    return ((1/Kt)*np.sqrt((PLsig1**2)+(PLsig1*E*((PLsig1/Hpr)**(1/npr))))-ProofS[2])
PLdelsig1 = optimize.brentq(f,ProofS[2], Kt*ProofS[2], args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
PLdelsig1

# Nonlinear root find for sigma of unload (factor of two expansion)
def f(PUsig1):
    return ((1/Kt)*np.sqrt(((PUsig1/2)**2)+((PUsig1*E)/2)*((PUsig1/(2*Hpr))**(1/npr)))-(ProofS[2]/2))
PUdelsig1 = optimize.brentq(f,ProofS[2], Kt*ProofS[2], args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
PUdelsig1

# Find sigma at new origin by subtracting unload stress from load stress
ProofSigma1 = PLdelsig1[0] - PUdelsig1[0]

# Find Epsilon at new origin with R.O. equations for loading and unloading
PLeps = (PLdelsig1[0]/E)+((PLdelsig1[0]/Hpr)**(1/npr))
PUeps = (PUdelsig1[0]/E) + 2*((PUdelsig1[0]/2/Hpr)**(1/npr))

# Subtract unload strain from load strain to find strain at new origin
ProofEpsilon1 = PLeps - PUeps

# New sigmax*eps_a values based on adjusting origin of S/S curve
SWTP1 = 1.852893797
SWTP2 = 0.168671136

# Nonlinear root find for Nf for each cycle
def f(NfP1):
    return (((sig_fp**2)/E)*((2*NfP1)**(2*b)) + (sig_fp*eps_fp*((2*NfP1)**(b+c))))-SWTP1
NfP1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)


def f(NfP2):
    return (((sig_fp**2)/E)*((2*NfP2)**(2*b)) + (sig_fp*eps_fp*((2*NfP2)**(b+c))))-SWTP2
NfP2 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

print(ProofEpsilon1)
print(ProofSigma1)
print(NfP1[0],NfP2[0])


# In[53]:


# Establish new Stress/Strain Origin after proof loading 130%

# Nonlinear root find for sigma of first load
def f(PLsig1):
    return ((1/Kt)*np.sqrt((PLsig1**2)+(PLsig1*E*((PLsig1/Hpr)**(1/npr))))-ProofS[3])
PLdelsig1 = optimize.brentq(f,ProofS[3], Kt*ProofS[3], args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
PLdelsig1

# Nonlinear root find for sigma of unload (factor of two expansion)
def f(PUsig1):
    return ((1/Kt)*np.sqrt(((PUsig1/2)**2)+((PUsig1*E)/2)*((PUsig1/(2*Hpr))**(1/npr)))-(ProofS[3]/2))
PUdelsig1 = optimize.brentq(f,ProofS[3], Kt*ProofS[3], args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
PUdelsig1

# Find sigma at new origin by subtracting unload stress from load stress
ProofSigma1 = PLdelsig1[0] - PUdelsig1[0]

# Find Epsilon at new origin with R.O. equations for loading and unloading
PLeps = (PLdelsig1[0]/E)+((PLdelsig1[0]/Hpr)**(1/npr))
PUeps = (PUdelsig1[0]/E) + 2*((PUdelsig1[0]/2/Hpr)**(1/npr))

# Subtract unload strain from load strain to find strain at new origin
ProofEpsilon1 = PLeps - PUeps

# New sigmax*eps_a values based on adjusting origin of S/S curve
SWTP1 = 1.663402872
SWTP2 = 0.064058188

# Nonlinear root find for Nf for each cycle
def f(NfP1):
    return (((sig_fp**2)/E)*((2*NfP1)**(2*b)) + (sig_fp*eps_fp*((2*NfP1)**(b+c))))-SWTP1
NfP1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)


def f(NfP2):
    return (((sig_fp**2)/E)*((2*NfP2)**(2*b)) + (sig_fp*eps_fp*((2*NfP2)**(b+c))))-SWTP2
NfP2 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

print(ProofEpsilon1)
print(ProofSigma1)
print(NfP1[0],NfP2[0])


# In[55]:


# Establish new Stress/Strain Origin after proof loading 140%

# Nonlinear root find for sigma of first load
def f(PLsig1):
    return ((1/Kt)*np.sqrt((PLsig1**2)+(PLsig1*E*((PLsig1/Hpr)**(1/npr))))-ProofS[4])
PLdelsig1 = optimize.brentq(f,ProofS[4], Kt*ProofS[4], args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
PLdelsig1

# Nonlinear root find for sigma of unload (factor of two expansion)
def f(PUsig1):
    return ((1/Kt)*np.sqrt(((PUsig1/2)**2)+((PUsig1*E)/2)*((PUsig1/(2*Hpr))**(1/npr)))-(ProofS[4]/2))
PUdelsig1 = optimize.brentq(f,ProofS[4], Kt*ProofS[4], args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
PUdelsig1

# Find sigma at new origin by subtracting unload stress from load stress
ProofSigma1 = PLdelsig1[0] - PUdelsig1[0]

# Find Epsilon at new origin with R.O. equations for loading and unloading
PLeps = (PLdelsig1[0]/E)+((PLdelsig1[0]/Hpr)**(1/npr))
PUeps = (PUdelsig1[0]/E) + 2*((PUdelsig1[0]/2/Hpr)**(1/npr))

# Subtract unload strain from load strain to find strain at new origin
ProofEpsilon1 = PLeps - PUeps

# New sigmax*eps_a values based on adjusting origin of S/S curve
SWTP1 = 1.534866653

# Nonlinear root find for Nf for each cycle
def f(NfP1):
    return (((sig_fp**2)/E)*((2*NfP1)**(2*b)) + (sig_fp*eps_fp*((2*NfP1)**(b+c))))-SWTP1
NfP1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)



print(ProofEpsilon1)
print(ProofSigma1)
print(NfP1[0])


# In[57]:


# Establish new Stress/Strain Origin after proof loading 150%

# Nonlinear root find for sigma of first load
def f(PLsig1):
    return ((1/Kt)*np.sqrt((PLsig1**2)+(PLsig1*E*((PLsig1/Hpr)**(1/npr))))-ProofS[5])
PLdelsig1 = optimize.brentq(f,ProofS[5], Kt*ProofS[5], args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
PLdelsig1

# Nonlinear root find for sigma of unload (factor of two expansion)
def f(PUsig1):
    return ((1/Kt)*np.sqrt(((PUsig1/2)**2)+((PUsig1*E)/2)*((PUsig1/(2*Hpr))**(1/npr)))-(ProofS[5]/2))
PUdelsig1 = optimize.brentq(f,ProofS[5], Kt*ProofS[5], args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)
PUdelsig1

# Find sigma at new origin by subtracting unload stress from load stress
ProofSigma1 = PLdelsig1[0] - PUdelsig1[0]

# Find Epsilon at new origin with R.O. equations for loading and unloading
PLeps = (PLdelsig1[0]/E)+((PLdelsig1[0]/Hpr)**(1/npr))
PUeps = (PUdelsig1[0]/E) + 2*((PUdelsig1[0]/2/Hpr)**(1/npr))

# Subtract unload strain from load strain to find strain at new origin
ProofEpsilon1 = PLeps - PUeps

# New sigmax*eps_a values based on adjusting origin of S/S curve
SWTP1 = 1.370500839

# Nonlinear root find for Nf for each cycle
def f(NfP1):
    return (((sig_fp**2)/E)*((2*NfP1)**(2*b)) + (sig_fp*eps_fp*((2*NfP1)**(b+c))))-SWTP1
NfP1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)



print(ProofEpsilon1)
print(ProofSigma1)
print(NfP1[0])


# In[65]:


# Nonlinear root find for Nf for each Cycle (-100%)

SWTP1 = 7.764084507
SWTP2 = 3.432083813
SWTP3 = 1.665889961
SWTP4 = 0.522363673

def f(NfP1):
    return (((sig_fp**2)/E)*((2*NfP1)**(2*b)) + (sig_fp*eps_fp*((2*NfP1)**(b+c))))-SWTP1
NfP1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)


def f(NfP2):
    return (((sig_fp**2)/E)*((2*NfP2)**(2*b)) + (sig_fp*eps_fp*((2*NfP2)**(b+c))))-SWTP2
NfP2 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

def f(NfP3):
    return (((sig_fp**2)/E)*((2*NfP3)**(2*b)) + (sig_fp*eps_fp*((2*NfP3)**(b+c))))-SWTP3
NfP3 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

def f(NfP4):
    return (((sig_fp**2)/E)*((2*NfP4)**(2*b)) + (sig_fp*eps_fp*((2*NfP4)**(b+c))))-SWTP4
NfP4 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

print(NfP1[0],NfP2[0],NfP3[0],NfP4[0])


# In[66]:


# Nonlinear root find for Nf for each Cycle (-110%)

SWTP1 = 8.056879359
SWTP2 = 3.593728137
SWTP3 = 1.773389156
SWTP4 = 0.576108721

def f(NfP1):
    return (((sig_fp**2)/E)*((2*NfP1)**(2*b)) + (sig_fp*eps_fp*((2*NfP1)**(b+c))))-SWTP1
NfP1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)


def f(NfP2):
    return (((sig_fp**2)/E)*((2*NfP2)**(2*b)) + (sig_fp*eps_fp*((2*NfP2)**(b+c))))-SWTP2
NfP2 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

def f(NfP3):
    return (((sig_fp**2)/E)*((2*NfP3)**(2*b)) + (sig_fp*eps_fp*((2*NfP3)**(b+c))))-SWTP3
NfP3 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

def f(NfP4):
    return (((sig_fp**2)/E)*((2*NfP4)**(2*b)) + (sig_fp*eps_fp*((2*NfP4)**(b+c))))-SWTP4
NfP4 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

print(NfP1[0],NfP2[0],NfP3[0],NfP4[0])


# In[67]:


# Nonlinear root find for Nf for each Cycle (-120%)

SWTP1 = 8.289705184
SWTP2 = 3.722265142
SWTP3 = 1.858870812
SWTP4 = 0.618845931

def f(NfP1):
    return (((sig_fp**2)/E)*((2*NfP1)**(2*b)) + (sig_fp*eps_fp*((2*NfP1)**(b+c))))-SWTP1
NfP1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)


def f(NfP2):
    return (((sig_fp**2)/E)*((2*NfP2)**(2*b)) + (sig_fp*eps_fp*((2*NfP2)**(b+c))))-SWTP2
NfP2 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

def f(NfP3):
    return (((sig_fp**2)/E)*((2*NfP3)**(2*b)) + (sig_fp*eps_fp*((2*NfP3)**(b+c))))-SWTP3
NfP3 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

def f(NfP4):
    return (((sig_fp**2)/E)*((2*NfP4)**(2*b)) + (sig_fp*eps_fp*((2*NfP4)**(b+c))))-SWTP4
NfP4 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

print(NfP1[0],NfP2[0],NfP3[0],NfP4[0])


# In[69]:


# Nonlinear root find for Nf for each Cycle (-130%)

SWTP1 = 8.479196109
SWTP2 = 3.82687809
SWTP3 = 1.928442123
SWTP4 = 0.653628643

def f(NfP1):
    return (((sig_fp**2)/E)*((2*NfP1)**(2*b)) + (sig_fp*eps_fp*((2*NfP1)**(b+c))))-SWTP1
NfP1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)


def f(NfP2):
    return (((sig_fp**2)/E)*((2*NfP2)**(2*b)) + (sig_fp*eps_fp*((2*NfP2)**(b+c))))-SWTP2
NfP2 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

def f(NfP3):
    return (((sig_fp**2)/E)*((2*NfP3)**(2*b)) + (sig_fp*eps_fp*((2*NfP3)**(b+c))))-SWTP3
NfP3 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

def f(NfP4):
    return (((sig_fp**2)/E)*((2*NfP4)**(2*b)) + (sig_fp*eps_fp*((2*NfP4)**(b+c))))-SWTP4
NfP4 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

print(NfP1[0],NfP2[0],NfP3[0],NfP4[0])


# In[70]:


# Nonlinear root find for Nf for each Cycle (-140%)

SWTP1 = 8.607732329
SWTP2 = 3.897839551
SWTP3 = 1.975634005
SWTP4 = 0.677222587

def f(NfP1):
    return (((sig_fp**2)/E)*((2*NfP1)**(2*b)) + (sig_fp*eps_fp*((2*NfP1)**(b+c))))-SWTP1
NfP1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)


def f(NfP2):
    return (((sig_fp**2)/E)*((2*NfP2)**(2*b)) + (sig_fp*eps_fp*((2*NfP2)**(b+c))))-SWTP2
NfP2 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

def f(NfP3):
    return (((sig_fp**2)/E)*((2*NfP3)**(2*b)) + (sig_fp*eps_fp*((2*NfP3)**(b+c))))-SWTP3
NfP3 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

def f(NfP4):
    return (((sig_fp**2)/E)*((2*NfP4)**(2*b)) + (sig_fp*eps_fp*((2*NfP4)**(b+c))))-SWTP4
NfP4 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

print(NfP1[0],NfP2[0],NfP3[0],NfP4[0])


# In[71]:


# Nonlinear root find for Nf for each Cycle (-150%)

SWTP1 = 8.772098143
SWTP2 = 3.988581587
SWTP3 = 2.035980669
SWTP4 = 0.707393365

def f(NfP1):
    return (((sig_fp**2)/E)*((2*NfP1)**(2*b)) + (sig_fp*eps_fp*((2*NfP1)**(b+c))))-SWTP1
NfP1 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)


def f(NfP2):
    return (((sig_fp**2)/E)*((2*NfP2)**(2*b)) + (sig_fp*eps_fp*((2*NfP2)**(b+c))))-SWTP2
NfP2 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

def f(NfP3):
    return (((sig_fp**2)/E)*((2*NfP3)**(2*b)) + (sig_fp*eps_fp*((2*NfP3)**(b+c))))-SWTP3
NfP3 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

def f(NfP4):
    return (((sig_fp**2)/E)*((2*NfP4)**(2*b)) + (sig_fp*eps_fp*((2*NfP4)**(b+c))))-SWTP4
NfP4 = optimize.brentq(f,1, 10**12, args=(), xtol=1e-10, maxiter=100, full_output=True, disp=True)

print(NfP1[0],NfP2[0],NfP3[0],NfP4[0])


# In[ ]:




