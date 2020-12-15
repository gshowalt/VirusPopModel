

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

from scipy.stats import ks_2samp


# define time, temperature scale
time = 5000
temp_list = np.linspace(-10, -1, 8)
t = np.linspace(1,time,1000)

# define runs, values for observed data
runs = 100

dfEric = pd.read_excel("ICEVBR_18Aug.xls")

Eric_VBR = dfEric.Virus/dfEric.Bacteria
RangelistB = []
RangelistM = []
RangelistD = []
RangelistP = []
RangelistBAvg = []
RangelistMAvg = []
RangelistDAvg = []
RangelistPAvg = []
simVBR  = []

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
    philo = np.random.uniform(-15,-5)
    phihi = np.random.uniform(philo,-5)



    # after we get a random range of values, we want to sample within that range of values to produce
    # a representative stretch of values to test if they actually reproduce the VBR.
    problem = {
        "num_vars" : 4,
        "names" : ['beta', 'mu', 'm', 'phi'],
        "bounds": [[betalo, betahi], 
                   [mulo, muhi], 
                   [deltalo, deltahi],
                   [philo, phihi]],
        "dists":['unif','unif','unif', 'unif']
    }

    param_values = saltelli.sample(problem,100,calc_second_order=True)
    

    # scale the parameters properly
    
    beta = param_values[:,0]
    mu = 10**param_values[:,1]
    delta = 10**param_values[:,2]
    phi = 10**param_values[:,3]

    # establish constant parameters
    alpha = (1.2e-7)*3**((-5-23)/10)#4.2e-7 at +8, or 1.2e-7 at lower temps, at -5 --> mu = 0.25/day = 0.01/hr = 1e-8
    #nutrient transfer coefficient to bacteria (ug/cell * hr)
    Q = 0.022
    #half saturation constant (ug/mL)
    d = 1e-8
    #constant of bacterial death (1/hr)
    
    # leak, lyse efficiencies here:
    g = 0.2
    n = 0.99
    
    
    # simulate VBRs
    simN = (n*1e-7*d*Q)/((alpha)*(g-1) + (n*1e-7)*(mu-d))
    simB = delta / (phi * (beta - 1))
    simV = ((mu*simN)/(simN + Q) - d)/(gamma * phi)
    
    simVBR = simV/simB
    
        # test the simulated VBRS against the real data
    result = ks_2samp(Eric_VBR, simVBR)
    if result[1] > 0.05:
        print('beta range is:', betalo, ',', betahi)        
        print('mu range is:', 10**(mulo), ',', 10**(muhi))
        print('delta range is:', 10**(deltalo), ',', 10**(deltahi), '\n')
        print('phi range is:', 10**(philo), ',', 10**(phihi), '\n')
        
        i += 1
        
        RangelistB.append([betalo, betahi])
        RangelistBAvg.append((betalo + betahi)/2)
        RangelistM.append([mulo, muhi])
        RangelistMAvg.append((mulo + muhi)/2)
        RangelistD.append([deltalo, deltahi])
        RangelistDAvg.append((deltalo + deltahi)/2)
        RangelistP.append([philo, phihi])
        RangelistPAvg.append((philo + phihi)/2)
    j += 1
    
print ('finished, total runs is:', j)
print ('ratio is:', i/j)


