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
from labellines import labelLine, labelLines
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
# define the function which includes the differential equations


def f2(s,t, temp, beta, mu, phi, delta, gamma):
    
    # first define the relative contact rate (RCR) and brine concentrating factor (BCF) by temp
    if temp < -1:
        RCR = 0.0716*temp**4 + 2.9311*temp**3 + 34.108*temp**2 + 45.826*temp + 3.5125 #Fit from Wells and Deming, 2006
        BCF = -0.0106 * temp **2 - 0.519 * temp + 0.2977
        sal = 32 * BCF
    else:
        RCR = 1
        sal = 32
  
    
    # scale adsorption rate by RCR to incorporate the sea ice 
    phi = phi * RCR    

    
    # SET PARAMETERS
    alpha = (1.2e-7)*3**((temp-23)/10)#4.2e-7 at +8, or 1.2e-7 at lower temps, at -5 --> mu = 0.25/day = 0.01/hr = 1e-8
    # alpha is a coefficient that we'd like to change with temperature? Or change eta?
    #nutrient transfer coefficient to bacteria (ug/cell * hr)
    Q = 0.022
    #half saturation constant (ug/mL)
    d = 1e-8
    #constant of bacterial death (1/hr)
    
    # leak, lyse efficiencies here:
    g = 0.1
    n = 0.99
    
    
    #gamma is a lysogeny value
    #gamma = 1 #-1/temp #*mu
    
    # set up solution matrix
    N = s[0]
    B = s[1]
    V = s[2]
  
    #systems of equations below
  
    dNdt = - alpha * (N / (N + Q)) * B + g * (alpha  * (N/(N+Q))*B) + (n * 1e-7 * (gamma) * phi * V * B) 
    if N < 0:
        N = 0
    dBdt = (mu) * (N/(Q + N)) * B - gamma * phi * V * B - d*B
    if B < 1:
        B = 1
    dVdt =  gamma*beta * B * phi*V - phi * V * B -  delta*V
    if V < 1:
        V = 1
   
    #print (mu, beta, phi, gamma)
    return [dNdt, dBdt, dVdt]
from scipy.stats import ks_2samp


# define runs, values for observed data
runs = 1

dfEric = pd.read_excel("ICEVBR_18Aug.xls")

Eric_VBR = dfEric.Virus/dfEric.Bacteria
RangelistB = []
RangelistM = []
RangelistD = []
RangelistBAvg = []
RangelistMAvg = []
RangelistDAvg = []

# in each run, we want to change the endpoints of the parameter distributions simultaneously
i = 0
j = 0



# the while loop is set to run until we find X sets of parameters that produce a distribution of 
# simulated VBRs which matches the distribution of real data with 95% acceptance rate

while i < runs:


    # define time, temperature scale
    time = np.random.uniform(5000,10000)
    temp_list = np.linspace(-12, -1, 8)
    t = np.linspace(1,time,2000)

    simVBR  = []



    
    
    # we want to include a step than randomly creates a range of values, since we're looking not for an indivudal 
    # value but instead looking for the right *range* of values.
    betalo = np.random.uniform(1,1000)
    betahi =  np.random.uniform(betalo,1000)
    mulo = np.random.uniform(-15,0)
    muhi = np.random.uniform(mulo,0)
    deltalo = np.random.uniform(-15,0)
    deltahi = np.random.uniform(deltalo,0)
    philo = np.random.uniform(-15,-1)
    phihi = np.random.uniform(philo,-1)
    #gammalo = np.random.uniform(0,1)
    #gammahi = np.random.uniform(gammalo,1)




    # after we get a random range of values, we want to sample within that range of values to produce
    # a representative stretch of values to test if they actually reproduce the VBR.
    problem = {
        "num_vars" : 2,
        "names" : ['mu', 'm'],
        "bounds": [[mulo, muhi], 
                   [deltalo, deltahi]],
        "dists":['unif','unif']
    }

    param_values = saltelli.sample(problem,100,calc_second_order=True)
    

    # scale the parameters properly
    
    beta = 100 #param_values[:,0]
    mu = 10**param_values[:,0]
    delta = 10**param_values[:,1]
    phi = 1e-10 #10**param_values[:,3]
    gamma = 1
    
    for temp in temp_list:
        if temp < -1:
            RCR = 0.0716*temp**4 + 2.9311*temp**3 + 34.108*temp**2 + 45.826*temp + 3.5125 #Fit from Wells and Deming, 2006
            BCF = -0.0106*temp**2 - 0.519*temp + 0.2977
            sal = 32 * BCF
        else:
            BCF = 1
            sal = 32

        nuts_init = 0.12 * BCF
        bact_init = 1e4 * BCF
        vir_init = 1e5 * BCF
        
        s0 = [nuts_init, bact_init, vir_init]
        for z in range(0,len(mu)):
            solx = odeint(f2, s0, t, args = (temp, beta, mu[z], phi, delta[z], gamma))

            nuts = solx[:,0]
            bact = solx[:,1]
            virus = solx[:,2]
            
            for x in bact:
                if x <= 0:
                    x = 5
            
            for x in virus:
                if x <= 0:
                    x = 5

            simVBR.append(virus/bact)    

    
    simVBR = np.concatenate(simVBR).ravel()

    # test the simulated VBRS against the real data
    result = ks_2samp(Eric_VBR, simVBR)
    if result[1] > 0.1:
        #print('beta range is:', betalo, ',', betahi)        
        print('mu range is:', 10**(mulo), ',', 10**(muhi))
        print('delta range is:', 10**(deltalo), ',', 10**(deltahi), '\n')
        
        i += 1
        
        RangelistB.append([betalo, betahi])
        RangelistBAvg.append((betalo + betahi)/2)
        RangelistM.append([mulo, muhi])
        RangelistMAvg.append((mulo + muhi)/2)
        RangelistD.append([deltalo, deltahi])
        RangelistDAvg.append((deltalo + deltahi)/2)

    j += 1
    
print ('finished, total runs is:', j)
print ('ratio is:', i/j)



"""# lets plot all the ranges as a set of ordered lines
#from itertools import izip
from scipy.stats import gaussian_kde


yax = []
yax1 = []
runs = i
    
for i in range(1,runs+1):
    yax.append([i,i])
    yax1.append(i)
    
sorted_lists = sorted(zip(RangelistB, RangelistM, RangelistD, RangelistP, RangelistG, RangelistMAvg, RangelistDAvg, RangelistPAvg, RangelistGAvg, RangelistBAvg ), reverse = True, key=lambda x: x[8])
RangelistB2, RangelistM2, RangelistD2, RangelistP2, RangelistG2, RangelistMAvg2, RangelistDAvg2, RangelistPAvg2, RangelistGAvg2, RangelistBAvg2 = [[x[i] for x in sorted_lists] for i in range(10)]

# And plot
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams.update({'font.size': 12.5})
#fig, axs = plt.subplots(1,2,figsize=(15,15))
fig = plt.figure(figsize=(20,10))
fig.tight_layout()
cm = plt.get_cmap('viridis')
#fig = plt.figure()
colorlist = [cm(1.*i/runs) for i in range(runs)]

ax = fig.add_subplot(251)
for i in range(runs):
    plt.plot(RangelistB2[i],yax[i], color = colorlist[i])
    plt.plot(RangelistBAvg2[i],yax1[i],color = colorlist[i], marker = 'o', markeredgecolor= 'k')
    plt.title('Burst Size')
    plt.xlabel('Burst Size')
    ax.set_xlim(0,1000)
    ax.set_yticklabels("")

ax1 = fig.add_subplot(252)
for i in range(runs):
    plt.plot(RangelistM2[i],yax[i], color = colorlist[i])
    plt.plot(RangelistMAvg2[i],yax1[i],color = colorlist[i], marker = 'o', markeredgecolor= 'k')
    plt.title('Growth Rate')
    plt.xlabel('Growth Rate Range (log 10)')
    ax1.set_yticklabels("")

ax2 = fig.add_subplot(253)
for i in range(runs):
    plt.plot(RangelistD2[i],yax[i], color = colorlist[i])
    plt.plot(RangelistDAvg2[i],yax1[i],color = colorlist[i], marker = 'o', markeredgecolor= 'k')
    plt.title('Decay Rate')
    plt.xlabel('Decay Rate Range (log 10)')
    ax2.set_yticklabels("") 
    
ax3 = fig.add_subplot(254)
for i in range(runs):
    plt.plot(RangelistP2[i],yax[i], color = colorlist[i])
    plt.plot(RangelistPAvg2[i],yax1[i],color = colorlist[i], marker = 'o', markeredgecolor= 'k')
    plt.title('Adsorption Rate')
    plt.xlabel('Adsorp. Rate Range (log 10)')
    ax3.set_yticklabels("")

ax4 = fig.add_subplot(255)
for i in range(runs):
    plt.plot(RangelistG2[i],yax[i], color = colorlist[i])
    plt.plot(RangelistGAvg2[i],yax1[i],color = colorlist[i], marker = 'o', markeredgecolor= 'k')
    plt.title('Lysis Rate')
    plt.xlabel('Lytic Fraction')
    ax3.set_yticklabels("")



ax5 = fig.add_subplot(256)
data = RangelistBAvg2
density = gaussian_kde(data)
xs = np.linspace(1,1000,1000)
density.covariance_factor = lambda : .25
density._compute_covariance()
plt.plot(xs,density(xs), 'k')
       
    
ax6 = fig.add_subplot(257)
data = RangelistMAvg2
density = gaussian_kde(data)
xs = np.linspace(-15,0,1000)
density.covariance_factor = lambda : .25
density._compute_covariance()
plt.plot(xs,density(xs), 'k')


ax7 = fig.add_subplot(258)
data = RangelistDAvg2
density = gaussian_kde(data)
xs = np.linspace(-15,0,1000)
density.covariance_factor = lambda : .25
density._compute_covariance()
plt.plot(xs,density(xs), 'k')

ax8 = fig.add_subplot(259)
data = RangelistPAvg2
density = gaussian_kde(data)
xs = np.linspace(-15,-5,1000)
density.covariance_factor = lambda : .25
density._compute_covariance()
plt.plot(xs,density(xs), 'k')

ax9 = fig.add_subplot(2,5,10)
data = RangelistGAvg2
density = gaussian_kde(data)
xs = np.linspace(0,1,1000)
density.covariance_factor = lambda : .25
density._compute_covariance()
plt.plot(xs,density(xs), 'k')
    
    
fig.tight_layout()
plt.subplots_adjust(top=0.85)

    
fig.suptitle('Full model parameter fitting', size = 30)"""
