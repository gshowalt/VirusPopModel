# This is us playing with simplified equations to determine the behavior of V*/B*

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

from scipy.stats import ks_2samp

runs = 50

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
    # we want to include a step than randomly creates a range of values, since we're looking not for an indivudal 
    # value but instead looking for the right *range* of values.
    betalo = np.random.uniform(1,1000)
    betahi =  np.random.uniform(betalo,1000)
    mulo = np.random.uniform(-15,0)
    muhi = np.random.uniform(mulo,0)
    deltalo = np.random.uniform(-15,0)
    deltahi = np.random.uniform(deltalo,0)




    # after we get a random range of values, we want to sample within that range of values to produce
    # a representative stretch of values to test if they actually reproduce the VBR.
    problem = {
        "num_vars" : 3,
        "names" : ['beta', 'mu', 'm'],
        "bounds": [[betalo, betahi], 
                   [mulo, muhi], 
                   [deltalo, deltahi]],
        "dists":['unif','unif','unif']
    }

    param_values = saltelli.sample(problem,100,calc_second_order=True)
    

    # scale the parameters properly
    
    beta = param_values[:,0]
    mu = 10**param_values[:,1]
    delta = 10**param_values[:,2]

    # simulate VBR
    simVBR = ((mu)*(beta - 1))/(delta)

    # test the simulated VBRS against the real data
    result = ks_2samp(Eric_VBR, simVBR)
    if result[1] > 0.05:
        print('beta range is:', betalo, ',', betahi)        
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



# lets plot all the ranges as a set of ordered lines
#from itertools import izip

yax = []
yax1 = []
    
for i in range(1,runs+1):
    yax.append([i,i])
    yax1.append(i)
    
sorted_lists = sorted(zip(RangelistB, RangelistM, RangelistD,RangelistMAvg, RangelistDAvg, RangelistBAvg ), reverse = True, key=lambda x: x[4])
RangelistB2, RangelistM2, RangelistD2, RangelistMAvg2, RangelistDAvg2,RangelistBAvg2 = [[x[i] for x in sorted_lists] for i in range(6)]


# And plot
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams.update({'font.size': 20})
#fig, axs = plt.subplots(1,2,figsize=(15,15))
fig = plt.figure(figsize=(20,5))
fig.tight_layout()
cm = plt.get_cmap('viridis')
#fig = plt.figure()
colorlist = [cm(1.*i/runs) for i in range(runs)]


ax = fig.add_subplot(131)
for i in range(runs):
    plt.plot(RangelistB2[i],yax[i], color = colorlist[i])
    plt.plot(RangelistBAvg2[i],yax1[i],color = colorlist[i], marker = 'o', markeredgecolor= 'k')
    plt.title('Burst Size')
    plt.xlabel('Burst Size')
    ax.set_yticklabels("")
    
ax1 = fig.add_subplot(132)
for i in range(runs):
    plt.plot(RangelistM2[i],yax[i], color = colorlist[i])
    plt.plot(RangelistMAvg2[i],yax1[i],color = colorlist[i], marker = 'o', markeredgecolor= 'k')
    plt.title('Growth Rate')
    plt.xlabel('Growth Rate Range (log 10)')
    ax1.set_yticklabels("")
    
ax3 = fig.add_subplot(133)
for i in range(runs):
    plt.plot(RangelistD2[i],yax[i], color = colorlist[i])
    plt.plot(RangelistDAvg2[i],yax1[i],color = colorlist[i], marker = 'o', markeredgecolor= 'k')
    plt.title('Decay Rate')
    plt.xlabel('Decay Rate Range (log 10)')
    ax3.set_yticklabels("")
    

fig.savefig('ParameterDistribution.png', dpi=300)


