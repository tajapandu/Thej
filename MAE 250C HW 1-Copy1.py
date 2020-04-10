#!/usr/bin/env python
# coding: utf-8

# In[171]:


import pylab
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import Series, DataFrame


# In[172]:


# 9.27
sigmaa = np.arange(10,3000,1)


# In[173]:


sigmafprime = 900
b = -0.102
Nf1 = (100+sigmaa)*sigmaa
Nf1= np.sqrt(Nf1)
Nf1 = Nf1/sigmafprime
Nf1 = Nf1**(1/b)
Nf1 = Nf1*.5


# In[174]:


Nf2 = sigmaa
Nf2 = Nf2/900
Nf2 = Nf2**(1/b)
Nf2 = Nf2*.5


# In[175]:


Nf3 = (sigmaa-100)*sigmaa
Nf3= np.sqrt(Nf3)
Nf3 = Nf3/sigmafprime
Nf3 = Nf3**(1/b)
Nf3 = Nf3*.5


# In[204]:


plt.plot(Nf1,sigmaa,'b',Nf2,sigmaa,'g',Nf3,sigmaa,'r')
plt.title('2024 T4 Al Stress Amplitude vs Cycle Life')
plt.xlabel(r'$N_f$')
plt.ylabel(r'$\sigma_a $'+ '(Mpa)')
plt.xscale('log')
plt.yscale('log')
plt.legend((r'$\sigma_m $'+'= 100 MPa',r'$\sigma_m $'+'= 0 MPa',r'$\sigma_m $'+'= -100 MPa'),
           loc='upper right', shadow=True)


# In[178]:


# 9.28
sigmam = np.arange(-400,800,1)
y1 = 1 - (sigmam/1172)
y2 = 1 - ((sigmam**2)/1373584)
y3 = 1 - (sigmam/1634)
y4 = 1 - (sigmam/1758)
datax = np.array([-207,-207,-207,207,207,207,414,414,414,414,414,414,621,623,621])
datay = np.array([1.122369,1.103515,1.179968,.910002,.913607,.958096,.672357,.785339,.745961,.763319,.846919,.865107,.613024,.635648,.68965])


# In[193]:


plt.plot(sigmam,y1,'r',sigmam,y2,'g',sigmam,y3,'b',sigmam,y4)
plt.scatter(datax,datay)
plt.title('Normalized Stress Amplitude-Mean Diagram')
plt.xlabel(r'$\sigma_m$')
plt.ylabel(r'$\frac{\sigma_a}{\sigma_ar}$')
plt.legend(('Goodman','Gerber','Morrow 1','Morrow 2'), loc='upper right', shadow=True)
plt.show()


# In[180]:


#9.29
sigmaar = np.array([948,834,703,631,579,524, 615.6298, 577.2954, 497.5661, 693.7773, 658.218, 585.4844, 511.7177, 511.7177, 
                    471.751, 647.2774, 577.2954, 507.0444, 471.2685, 436.3943, 365.1137])
Nf = np.array([222, 992, 6004, 14130, 43860, 132150, 73780, 83810, 567590, 31280, 50490, 84420, 431170, 730570, 445020, 
               45490, 109680, 510250, 208030, 193220, 901430 ])
Nfmodel = np.arange(10,950000,10)
sigmaarmodel = 1758*((2*Nfmodel)**-0.0977)


# In[202]:


plt.plot(Nfmodel,sigmaarmodel,'g')
plt.scatter(Nf, sigmaar)
plt.title('SWT ' + r'$\sigma_{ar}$ ' + 'vs ' + r'$N_f $')
plt.ylabel(r'$\sigma_{ar}$')
plt.xlabel(r'$N_f$')
plt.xscale('log')
plt.yscale('log')


# In[203]:


sigmaw = np.array([948, 834,703,631,579,524,532.2498,494.6801,416.9349083,635.1864739,599.8504515, 527.6680997,454.6393618, 
                   454.6393618, 417.1512979, 617.0850459, 547.2210773, 477.1251196, 503.1033092, 468.270807, 
                   397.0853035])
plt.plot(Nfmodel,sigmaarmodel,'r')
plt.title('Walker ' + r'$\sigma_{ar}$ ' + 'vs ' + r'$N_f $')
plt.scatter(Nf, sigmaw)
plt.xlabel(r'$N_f$')
plt.ylabel(r'$\sigma_{ar}$')
plt.xscale('log')
plt.yscale('log')


# In[ ]:





# In[ ]:




