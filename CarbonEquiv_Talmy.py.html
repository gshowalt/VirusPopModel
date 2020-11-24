#!/usr/bin/env python
# coding: utf-8


# Import all modules
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import matplotlib.tri as tri
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
from matplotlib.ticker import LogFormatter 
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


### ---- define the equation with system of ODEs ---- ###
def f2(s,t, leak, lyse, temp):
    
    if temp < -1:
        RCR = 0.0716*temp**4 + 2.9311*temp**3 + 34.108*temp**2 + 45.826*temp + 3.5125 #Fit from Wells and Deming, 2006
        BCF = -0.0106 * temp **2 - 0.519 * temp + 0.2977
        sal = 32 * BCF
    else:
        RCR = 1
        sal = 32
  
    #build burst size curve as fxn of Temp
    mux = 1e-2
    betx = 1
    phix = 1e-5
    gamx = 1
    
    
    beta = betx*(0.0064 * temp**3 - 0.3047 * temp ** 2 + 0.7701 * temp + 93.605)
    """ beta_std = 0.0095 * temp **3 - 0.5184 * temp**2 + 2.2456 * temp + 126.59
    if beta_std < 0:
        beta_std = 0.
    beta = np.random.normal(beta_mu, beta_std)"""

    #build growth curve as fxn of Temp
    #mu = mux*(2e-5*temp**3 + 0.0008 * temp **2 + 0.0091 * temp + 0.0386)
    #mu = 3e-6*temp**4 + 0.0001*temp**3+0.0014*temp**2 + 0.0092 * temp +0.0333
    mu = 0.0441*np.exp(0.4991*temp)
    """mu_std = 0.1*2e-5*temp**3 + 0.0009 * temp **2 + 0.0144 * temp + 0.0818
    if mu_std<0:
        mu_std = 0.001
    mu = np.random.normal(mu_mu, mu_std)"""

    #build adsorption curve as fxn of Salinity (which is function of temp)
    #phi = phix * -1e-11*sal**2 +4e-9*sal - 9e-8
    phi = phix * (6e-13 * temp **5 - 2e-11 * temp ** 4 + 1e-10 * temp ** 3 + 3e-9 * temp ** 2 - 3e-8 * temp + 5e-8)
    """phi_std = -2e-11*sal**2 + 4e-9*sal - 9e-8
    if phi_std < 0:
        phi_std = 0
    phi = np.random.normal(phi_mu, phi_std)"""
    
    if mu <= 0:
        mu = 0.000
    if beta < 0:
        beta = 1
    if phi < 0:
        phi = 1e-15
    
    
    
    phi = phi * RCR
    #print (phi)
    

    
    #print ((phi))
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
    N = s[0]
    B = s[1]
    V = s[2]
    P = s[3]
    #systems of equations below
    if N < 0:
        N = 0
    if B < 1:
        B = 1
    if V < 1:
        V = 1
    dNdt = - alpha * (N / (N + Q)) * B + g * (alpha  * (N/(N+Q))*B) + (n * 1e-7 * (gamma) * phi * V * B)
    if N < 0:
        N = 0
    #nutrient term
    dBdt = (mu) * (N/(Q + N)) * B - gamma * phi * V * B - d*B
    if B < 1:
        B = 1
    dVdt =  gamma*beta * B * phi*V - phi * V * B -  m*V
    if V < 1:
        V = 1
    #virus term
    #dPdt = (g * (0.0083*1e-7))*B + (n * 1e-7 * phi * V * B*RCR) + 1e-10*m*V + 1.0e-7*d*B - (P/(P+Q))*alpha * B
    dPdt = g * alpha  * (N/ (N+Q))*B + n * 1e-7 * (gamma)*phi*B*V
    #POM term
   
    # according to Jover, 2014 - virus has 0.02 to 0.05 fg carbon/virion => translate into ug Carbon = 5e-11
    VCarbonEQ = 5e-11
    BCarbonEQ = 1e-7 #from Bionumbers
    
    
    rv = 90 #virus radius
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


### ---- Set up condition and run the equation ---- ##
time = 5000
temp_list = [-12.5,-10, -8, -6, -4, -2]
t = np.linspace(1,time,1000)


DOMX = []
DOMA = []
DOMB = []
DOMC = []
DOM1 = []
DOM10 = []
DOM100 = []

vari = [0, 0.001, 0.01, 0.10, 0.50, 0.80, 1]

xend = []
aend = []
bend = []
cend = []
dend = []
eend = []
fend = []


#vari = np.arange(0,1000,0.001)
RCRlist = []
Mulist = []

vari2 = np.arange(0,1.1,0.1)
vari1 = np.arange(0,0.4,0.1)
count = 0
plt.rcParams["font.family"] = "sans-serif"
fig2 = plt.figure(figsize=(30,15))
fig2.tight_layout()
plt.rcParams.update({'font.size': 20})
endvals1 = []
endvals2 = []
endvals3 = []
endvals4 = []
Burstlist = []
Adsorplist = []

for xx in temp_list:
    temp = xx
    count +=1
    xend = []
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
        
    #print ("Temp is:", temp, "and RCR is:", RCR)
    
    for i in vari1:
        for j in vari2:
            
            s0=[0.12*BCF,1e4*BCF, 1e5*BCF,0,0,0,0,0]
            s = odeint(f2,s0,t, args = (i,j, temp))
            xend.append(sum(s[:,3]))

    #print (xend)     
    length =  len(xend)
    #print (length)
    y = list(islice(cycle(vari1),length))
    vari_axis = []

    for i in vari1:
        vari_time = [i]*11
        vari_axis  = vari_axis + vari_time

    y1 = s[:,4]/(0.12*BCF)
    y2 = s[:,5]/(0.12*BCF)
    y3 = s[:,6]/(0.12*BCF)
    y4 = s[:,7]/(0.12*BCF)
    
   
    plt.subplot(3, 3, count)

    
    colors1 = ['cadetblue', '#FF6F61'] #, 'darkblue']
    plt.stackplot(t,y2,y3, colors = colors1,labels=['To Virus','To Lysate'])
    plt.legend(loc='lower right')

    plt.xlabel('Temperature: {} (˚C)'.format(temp))
    plt.yscale('log')
    plt.ylabel('% Initial Nutrient')


    
    
    endvals1.append(y1[-1])
    endvals2.append(y2[-1])
    endvals3.append(y3[-1])
    endvals4.append(y4[-1])
    RCRlist.append(RCR)
    Mulist.append(mu)
    beta = 1*(0.0064 * temp**3 - 0.3047 * temp ** 2 + 0.7701 * temp + 93.605)
    Burstlist.append(beta)
    phi = RCR* 1 * (6e-13 * temp **5 - 2e-11 * temp ** 4 + 1e-10 * temp ** 3 + 3e-9 * temp ** 2 - 3e-8 * temp + 5e-8)
    Adsorplist.append(phi)

#fig1.text(0.5, 0.36, 'fraction of exudate to DOC pool',ha='center',fontsize=10)
#fig1.text(0.08, 0.65, 'fraction of lytic material to DOC pool', va='center', rotation='vertical',fontsize=10)



plt.subplots_adjust(hspace = 0.25)


"""plt.subplot(3,1,3)
slope = [-1.6,-0.14,-0.15,-0.22,-0.6,-2.15]
#note that slopes were calculated using web plot digitizer and data is available on my desktop
plt.plot(temp_list, slope, 'k')
plt.title('Relative Lytic Control', fontsize = 15)
plt.ylabel('Contour Slope (0 = Max Lytic)', fontsize=10)
plt.xlabel('Temperature (˚C)', fontsize=10)"""
#ax1.text(textx,texty,"g = 0", color = "white")
#ax1.set_title("DOM regenerated")
#ax1.set_ylim([0.01,2])
#ax1.set_yscale('log')

fig1.suptitle("Cumulative organic carbon recycled as a function of driving source ",fontsize=20)

#plt.stackplot(temp,endvals2,endvals3,labels=['To Virus','To Lysate'])




### ---- plot as a function of temp rather than time ---- ##
plt.rcParams["font.family"] = "sans-serif"
fig2 = plt.figure(figsize=(20,10))
fig2.tight_layout()

plt.rcParams.update({'font.size': 20})
ax1 = plt.stackplot(temp_list, endvals2, endvals3, colors = colors1) #, labels=['To Virus','To Lysate', 'Cell exudate'])
#ax1 = plt.plot(temp_list, Burstlist)


#plt.legend(loc='upper left')
plt.ylabel('Carbon Movement')
plt.xlabel('Temperature')


# In[105]:


#fig1.savefig('CE_Grid.jpeg',  bbox_inches="tight", dpi=300,transparent=True)
fig2.savefig('CE_GrowthratevTemp.jpeg',  bbox_inches="tight", dpi=300,transparent=True)


# In[ ]:




