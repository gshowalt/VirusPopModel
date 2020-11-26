
# importing all modules
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
from matplotlib.ticker import LogFormatter 

from collections import Counter

from functools import wraps

import csv
import sys

import itertools
from itertools import islice, cycle, chain

import scipy as sp
from scipy.interpolate import griddata
from scipy import interpolate
from scipy.integrate import odeint
from scipy.stats import pareto
from scipy.stats import loguniform

import seaborn as sns
import pandas as pd

import statistics as stats
import lhsmdu

from math import nan

from SALib.sample import saltelli, latin, ff
from SALib.analyze import sobol

import random


# In[127]:


# define the function which includes the differential equations
# this was adapted from the leak/lyse experiment so I just left that in and set it to a final value later

def f2(s,t, leak, lyse, temp):
    
    # first define the relative contact rate (RCR) and brine concentrating factor (BCF) by temp
    if temp < -1:
        RCR = 0.0716*temp**4 + 2.9311*temp**3 + 34.108*temp**2 + 45.826*temp + 3.5125 #Fit from Wells and Deming, 2006
        BCF = -0.0106 * temp **2 - 0.519 * temp + 0.2977
        sal = 32 * BCF
    else:
        RCR = 1
        sal = 32
  
    # these are our scaling factors for the temperature-dependent parameter distributions
    mux = 1       # for growth rate
    betx = 1      # for burst size
    phix = 1e-5   # for adsorption rate
    gamx = 1      # for lytic fraction
    
    # Temp-dependent parameter distribution for burst size
    beta = betx*(0.0064 * temp**3 - 0.3047 * temp ** 2 + 0.7701 * temp + 93.605)
    # also parameterized as a curve with a standard deviation (std) for other experiments
    # but here was simply a set curve for reproducibility
    """ beta_std = 0.0095 * temp **3 - 0.5184 * temp**2 + 2.2456 * temp + 126.59
    if beta_std < 0:
        beta_std = 0.
    beta = np.random.normal(beta_mu, beta_std)"""

    # Temp-dependent parameter distribution for growth rate 
    # (we had two different distributions, but I went with the exponential one)
    # mu = mux*(2e-5*temp**3 + 0.0008 * temp **2 + 0.0091 * temp + 0.0386)
    # mu = 3e-6*temp**4 + 0.0001*temp**3+0.0014*temp**2 + 0.0092 * temp +0.0333
    mu = 0.0441*np.exp(0.4991*temp) 
    """mu_std = 0.1*2e-5*temp**3 + 0.0009 * temp **2 + 0.0144 * temp + 0.0818
    if mu_std<0:
        mu_std = 0.001
    mu = np.random.normal(mu_mu, mu_std)"""

    # Temp-dependent parameter distribution for adsorption rate 
    # I also tried it as a function of salinity (immediately below), but chose temp for consistency
    #phi = phix * -1e-11*sal**2 +4e-9*sal - 9e-8
    phi = phix * (6e-13 * temp **5 - 2e-11 * temp ** 4 + 1e-10 * temp ** 3 + 3e-9 * temp ** 2 - 3e-8 * temp + 5e-8)
    """phi_std = -2e-11*sal**2 + 4e-9*sal - 9e-8
    if phi_std < 0:
        phi_std = 0
    phi = np.random.normal(phi_mu, phi_std)"""
    
    # set conditions for when curve goes below zero
    if mu <= 0:
        mu = 0.000
    if beta < 0:
        beta = 1
    if phi < 0:
        phi = 1e-15
    
    # now we want to scale adsorption rate by RCR to incorporate the sea ice 
    phi = phi * RCR    

    
    # SET PARAMETERS
    alpha = 1.2e-7*3**((temp-23)/10)#4.2e-7 at +8, or 1.2e-7 at lower temps, at -5 --> mu = 0.25/day = 0.01/hr = 1e-8
    # alpha is a coefficient that we'd like to change with temperature? Or change eta?
    #nutrient transfer coefficient to bacteria (ug/cell * hr)
    Q = 0.022
    #half saturation constant (ug/mL)
    d = 1e-8
    #constant of bacterial death (1/hr)
    m = 1e-6
    #constant of viral decay (1/hr)
    g = leak
    #POM transfer coefficient from bacteria (ug/cell*hr)
    n = lyse
    #POM transfer coefficient from viral lysis ug/[burst]cell
    #gamma is a lysogeny value
    gamma = 1 #-1/temp #*mu
    
    # set up solution matrix
    N = s[0]
    B = s[1]
    V = s[2]
    P = s[3]
    
    #systems of equations below
  
    dNdt = - alpha * (N / (N + Q)) * B + g * (alpha  * (N/(N+Q))*B) + (n * 1e-7 * (gamma) * phi * V * B)
    if N < 0:
        N = 0
    dBdt = (mu) * (N/(Q + N)) * B - gamma * phi * V * B - d*B
    if B < 1:
        B = 1
    dVdt =  gamma*beta * B * phi*V - phi * V * B -  m*V
    if V < 1:
        V = 1
    #dPdt = (g * (0.0083*1e-7))*B + (n * 1e-7 * phi * V * B*RCR) + 1e-10*m*V + 1.0e-7*d*B - (P/(P+Q))*alpha * B
    dPdt = g * alpha  * (N/ (N+Q))*B + n * 1e-7 * (gamma)*phi*B*V
   
    # according to Jover, 2014 - virus has 0.02 to 0.05 fg carbon/virion => translate into ug Carbon = 5e-11
    VCarbonEQ = 5e-11
    BCarbonEQ = 1e-7 #from Bionumbers
    
    # building the carbon equivalent for viruses, lysate as per Talmy et al 2019
    rv = 90 #virus radius (nm)
    Qv = (41 * (rv - 2.5)**3 + 130*(7.5*(rv)**2 - 18.74 * rv + 15.63)) * (10e6/(6.022 * 10**23)) # virus carbon eq
    phiEQ = (phi)/(Qv) 
    Qh =  1e-7
    etav = beta * (Qv/Qh)
    
    TotalVCarbon = (phiEQ * (gamma) * (V*VCarbonEQ) * (B*BCarbonEQ))
    VirusCarbon = etav * (phiEQ * (gamma) * (V*VCarbonEQ) * (B*BCarbonEQ))
    LysateCarbon = (1-etav)*(phiEQ * (gamma) * (V*VCarbonEQ) * (B*BCarbonEQ))
    LeakCarbon = g * (alpha  * (N/(N+Q))*B)

    
    #print (mu, beta, phi, gamma)
    return [dNdt, dBdt, dVdt, dPdt, TotalVCarbon, VirusCarbon, LysateCarbon, LeakCarbon]


# In[133]:


# define time, temperature scale
time = 5000
temp_list = [-12.5,-10, -8, -6, -4, -2]
t = np.linspace(1,time,1000)

# set up empty matricies
DOMX = []
DOMA = []
DOMB = []
DOMC = []
DOM1 = []
DOM10 = []
DOM100 = []

RCRlist = []
Mulist = []
endvals1 = []
endvals2 = []
endvals3 = []
endvals4 = []
Burstlist = []
Adsorplist = []

count = 0
plt.rcParams["font.family"] = "sans-serif"
fig1 = plt.figure(figsize=(20,15))
fig1.tight_layout()
plt.rcParams.update({'font.size': 15})

for xx in temp_list:
    temp = xx
    count +=1
    mu = 0.0441*np.exp(0.4991*temp)
    gamma = 1
    #print ("gamma is:", gamma, "and mu is:", mu)
    if temp < -1:
        RCR = 0.0716*temp**4 + 2.9311*temp**3 + 34.108*temp**2 + 45.826*temp + 3.5125 #Fit from Wells and Deming, 2006
        BCF = -0.0106 * temp **2 - 0.519 * temp + 0.2977
        sal = 32 * BCF
    else:
        BCF = 1
        sal = 32
        
    s0=[0.12*BCF,1e4*BCF, 1e5*BCF,0,0,0,0,0]
    s = odeint(f2,s0,t, args = (0.4,0.99, temp))
    xend.append(sum(s[:,3]))
 
    
    y1 = s[:,4]/(0.12)
    y2 = s[:,5]/(0.12)
    y3 = s[:,6]/(0.12)
    y4 = s[:,7]/(0.12)
    
   
    plt.subplot(3, 3, count)

    
    colors1 = ['cadetblue', '#FF6F61'] #, 'darkblue']
    plt.stackplot(t,y2,y3, colors = colors1,labels=['To Virus','To Lysate'])
    plt.legend(loc='lower right')

    plt.xlabel('Temperature: {} (˚C)'.format(temp))
    plt.yscale('log')
    plt.ylabel('% Initial Nutrient')


    
    # take last value of each returned number for the temp-dependent plot 
    endvals1.append(y1[-1])
    endvals2.append(y2[-1])
    endvals3.append(y3[-1])
    endvals4.append(y4[-1])
    
    # make lists of calculated temp-dependent parameters if we want to plot against them alter
    RCRlist.append(RCR)
    Mulist.append(mu)
    beta = 1*(0.0064 * temp**3 - 0.3047 * temp ** 2 + 0.7701 * temp + 93.605)
    Burstlist.append(beta)
    phi = RCR* 1 * (6e-13 * temp **5 - 2e-11 * temp ** 4 + 1e-10 * temp ** 3 + 3e-9 * temp ** 2 - 3e-8 * temp + 5e-8)
    Adsorplist.append(phi)



plt.subplots_adjust(hspace = 1)
fig1.suptitle("Cumulative organic carbon recycled into Virions or Lysate ",fontsize=15)

# Plot as a funciton of temperature
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams.update({'font.size': 20})
fig2 = plt.figure(figsize=(10,5))
fig2.tight_layout()


endvals1_b = [i/max(endvals1) for i in endvals1]
endvals2_b = [i/max(endvals2) for i in endvals2]
endvals3_b = [i/max(endvals3) for i in endvals3]
endvals4_b = [i/max(endvals4) for i in endvals4]

#ax1 = plt.stackplot(temp_list, endvals2_b, endvals3, colors = colors1) #, labels=['To Virus','To Lysate', 'Cell exudate'])
#ax1 = plt.plot(temp_list, Burstlist)
plt.plot(temp_list,endvals2_b, c = 'cadetblue', marker = 'o', markeredgecolor='white', markersize=15, label='to Virions')
plt.plot(temp_list, endvals3_b, c = '#FA7268', marker = 'o', markeredgecolor='white', markersize=15, label='to Lysate') 

plt.xlabel('Temperature (˚C)')
plt.ylabel('Carbon Flow (Relative to Maximum)')
plt.legend(loc='lower right')
fig2.suptitle("Cumulative organic carbon recycled into \nVirions or Lysate as a function of temperature\n",fontsize=15)




# In[88]:
#fig1.savefig('CE_Grid_withRCR_runaway.jpeg',  bbox_inches="tight", dpi=300,transparent=True)
#fig2.savefig('CE_Temp_noRCR_line.jpeg',  bbox_inches="tight", dpi=300,transparent=True)



