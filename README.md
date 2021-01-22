# Modeled virus-bacteria dynamics in sea ice brines

### Abstract 

Observations indicating high concentrations of viruses within sea-ice brines relative to bacteria indicate that bacteriophage may play an outsized role in shaping host community structure and recycling organic matter within sea ice. However, while these high virus-to-bacteria ratios may be a function of active viral production and release, they may also be the result of slow virion decay due to the low temperature, high salinity conditions of sea-ice brines. Here, we present a mathematical model of virus–host interactions within sea-ice brines, one that uses viral and bacterial abundance data to interpret which mechanisms may most influence population dynamics. We also assess the existence and potential impact of a viral shunt within sea ice. Data from both field samples and laboratory isolates were used to achieve most likely parameter distributions for *in situ* communities, constraining the model to observed dynamics. Deeper understanding of the viral impact on sea-ice communities and the recycling of dissolved organic carbon within the ice will help to understand how these largely heterotrophic communities maintain activity under extreme conditions and further constrain Arctic carbon budgets through the seasons.

# Introduction

Bacteriophage, viruses which infect bacteria, are numerically dominant within the global oceans and contribute 20-40% of daily bacterial mortality in pelagic environments (Suttle, 2007). As a result, viruses play a key role in microbial population dynamics and the movemment or organic matter within the global oceans. Notably, viral lysis of prokayotes releases organic matter, short-circuting the microbial loop in a process known as the viral shunt. The viral shunt is suggested to recycle as much as 25% of all carbon fixed by autotrophs in the pelagic ocean, facilitating heterotrophic microbial growth and activity (Wilhelm & Suttle, 1999; Weitz & Wilhelm, 2012).

When seawater freezes, both viruses and bacteria are incorporated into a porous matrix of briney liquid. Within this brine network of sea ice, microbial communities experience  extremes of temperature (from near –40˚C to 0˚C) and salinity (from near 240 ppt to 0 ppt salts), as well as physical concentration leading to cell densities up to 10^7 cells per milliliter of brine (Ewert & Deming, 2013; Boetius et al., 2015). Field observations of viral and bacterial abundances in the form of virus-to-bacteria ratios within sea ice, which approach a ratio of 10,000-to-1, can greatly exceed those of source seawater (Collins & Deming, 2011; Figure 1), which often range between 1 and 100 viruses per bacterium (Knowles et al., 2016; Wigington, et al., 2016; Parikka et al., 2017). These very high ratios have in turn led to the suggestion that viruses must be produced within sea-ice brines at rates well beyond those observed in seawater (Collins & Deming, 2011; Wigington et al., 2016). The primary argument supporting this suggestion is that tight physical coupling between virus and host within a brine pore favors increased infection rates. Indeed, as a result of the physical freeze-concentration effect and reduced viscosity of subzero brines, the relative contact rate of viruses and bacteria within brines has been calculated to be up to 1,000 times greater than in seawater, depending on temperature (Wells & Deming, 2006a). 

![Fig](https://github.com/gshowalt/VirusPopModel/blob/main/VBRplot_2.png)

*Figure 1: Density plot of observed virus-to-bacteria ratios in sea ice brines (reds) and seawater (blues).*


**We hypothesize that, as a result of physical concentration enhancing viral infection, viruses within sea ice brines may play an out-sized role in facilitating carbon cycling by contributing dissolved organic carbon to sea ice brine DOC pools in the form of lytic material.**

Here, we have built a simple population dynamic model using ordinary differential equations to demonstrate potential and constrain rates of viral infection and virally-mediated carbon cycling within sea ice brines. 


# Methods

### The Equations ###

Equations in this system were modeled on those from given in [*Quantitative Viral Ecology* by Joshua Weitz (2015)](https://press.princeton.edu/books/hardcover/9780691161549/quantitative-viral-ecology) These modles use simple differential equations to build interacting populations of bacteria (B) and viruses (V) by modeling their change as a function of time and nutrient concentration (N), such as:

Bacterial Population
> dBdt = growth - infection - death 

> dBdt = µ * N -  φ * V * B - δb * B 

Virus Population
> dVdt = production (*lysis*) - adsorption - decay 

> dVdt = β * φ * V * B -  φ * V * B - δv * V 


Which were then coded into Python 3:
```
dNdt = (-alpha * (N / (N + Q)) * B) + (g * (alpha  * (N/(N+Q))*B)) + (n * 1e-7 * phi * V * B)
dBdt = ((mu) * (N/(Q + N)) * B) - (phi * V * B) - d*B
dVdt =  (beta * B * phi * V) - (phi * V * B) -  (m * V)
```

### Assumptions ###
These equations differ from those by Weitz or similar open-water models in that they _do not contain any nutrient inflow_ - instead, they consider the sea ice brine pore as a closed system. As a result - in the equations given above - nutrients are consumed in the given equations and only recycled through exudate and lysis (not through bacterial death).

We made use of an existing parameterize of relative contact rate (RCR) produced by Wells and Deming (2006b), which calcualted how often viruses should encounter bacteria within sea ice as a result of increased physical concentration and decreased particle diffusion compared to warmer seawater conditions. This figure is given in file [WellsRCR](https://github.com/gshowalt/VirusPopModel/blob/main/WellsRCR.png).

Bacterial Population
> dBdt = growth - infection - death 

> dBdt = µ * N -  **RCR** * φ * V * B - δb * B 

Virus Population
> dVdt = production (*lysis*) - adsorption - decay 

> dVdt = β * **RCR** * φ * V * B -  **RCR** * φ * V * B - δv * V 


# Results & Discussion

## 1. Basic system behavior

When the system is run with values taken from literature, given in [Table 1](https://github.com/gshowalt/VirusPopModel/blob/main/Table1.png), we can see time-dependent dynamics. 

!(Fig)[https://github.com/gshowalt/VirusPopModel/blob/main/PhasePortrait_noannot_21Jan.png

These relationships can be easily visualized when looking at time-dependant behavior of the system between 0 and 10,000 hours in the phase plane, as plotted below.

First, this figure represents the phase portrait of bacterial and viral populations across the temperature spectrum for the equations run without relative contact rate coefficient (left panel, plasma colorscheme) and with relative contact rate coefficient (right panel, viridis colorscheme). Phase plots begin at the un-marked end of the line and run down the plot, ending at the solid black point (10,000 hours). Faded colored plots represent the same system without lytic recycling, while the gray spotted line shows VBR.



The pattern of these plots generally reflects the expected semi-cyclic behavior of both bacterial and viral populations within the system: a increase in bacterial population (movement -->) leads to a growth in viral population (movement ^), which inturn leades to a crash of the bacterial population (<--), creating a spiral/nautilis-like pattern. Two major execptions are obvious, both occuring at low temperature. In both plots, there is little change in bacterial population at low temperatures due to growth, indicated by little horizontal movement of the plots (especially towards positive bacterial growth). More noticably, without physical concentration increasing viral-host encounters (left panel), there is minimal infection (no vertical change) and indeed little change in either bacterial or viral populations below –6˚C.This is not observed in the right panel, where physical concentration allows for enhanced bacterial infection and therefor higher viral populations. Comparing these figures to the original VBR density plot (above), we see that including RCR creates time-dependent populations which largely overlap the observed VBR space.

## 2. Organic material cycling within the system

When we compared bacterial populations to nutrient/organic material concentration in the phase plane, we can again see the importance of physical concentration, as well as the effects of lytic recycling within the system (below). Again, the left panel represents the system run without physical concentration due to freezing, while the right panel represents the system with physical concentration. In both panels, the faded colored lines represent the same systems without lytic recycling.

In the left panel with no physical concentration, the lack of viral infection is obvious at temperatures below –6˚C. Minimal horiztonal change in these lines indicates little growth or loss of cells due to death or lysis. Conversely, we can see the signal of viral infection in the right panel with physical concentration coefficient included: as lines move from t0 toward tend, the first fall vertically (representing initial uptake) with a horizontal component (representing growth), but sharply move left (bacterial loss) with a vertical component (lysis). The difference between the solid and spotted lines indicate the relative input of lytic recycling compared to other sources of organic recycling (viral decay, bacteria exudation). It is clear when comparing final nutrient concentrations with and without lytic recycling (on solid or faded dotted lines, respetively) that the nutrient concentration after 10,000 hours is higher with lytic recycling included, and differentially so comparing the left and right panels. 

![Phase Portrait, Nuts](https://github.com/gshowalt/VirusPopModel/blob/main/PhasePortrait_nuts_17Jan.png)

We can pull out the relative contribution of lytic recycling by comparing the final nutrient concentration the runs at each temperature, subtracting the nutrient concentration of the run _without lytic recycling_ to see what percent of the total recycled nutrients comes from lytic material after 10,000 hours, as shown in the figure below. We can see that lysis contributes a large fraction of the total recycled material for both systems (with and without physical concentration), but the relative importance of lytic material depends on temperature. At colder temperatures (below –9˚C), the system with physical concentration factor shows a higher relative input of lytic material, which is consistent with our observations that there is less bacterial infection at lower temperatures if no physical concentration is included. This system is modeled assuming nearly all lytic material (99%) is recycled, and cells exude 10% of all material they take up back into a pool of recycled organic material.

![Phase Portrait, Nuts](https://github.com/gshowalt/VirusPopModel/blob/main/Lytic%20contribution_17Jan.png)

If we were to look at relative lytic contribution to recycled material across all possible parameter space, we can co-vary lytic and exudate material fractions to invesitgate if this patterns holds. Indeed, we see that at the "moderate" temperature ranges ( between –12˚C and –4˚C), lysis is more important relative to exudate, while the converse is true at endpoint temperatures (-13˚C and -2˚C).

** INSERT FIGURE FROM DISSERTATION ***

With this in mind, we can translate our bacterial and viral populations into carbon equivalents to estimate relative and absolute amount of organic carbon recycled.

The code for this step is uploaded to the repository as [CarbonEquiv_noRCR.py](https://github.com/gshowalt/VirusPopModel/edit/main/Code/CarbonEquiv_Talmy.py) .

Here's our figure of carbon cycling (1 = maximum carbon) as a function of temperature **WITHOUT** Brine Concentrating Factor

![Practice Text for Sizing](https://github.com/gshowalt/VirusPopModel/blob/main/CE_Temp_noRCR_line.jpeg)

Secondly, we apply the **physical concentration parameter** to the model. This physical concentration, a result of tightly constricted pore space within sea ice crystals, has been suggested to increase the _relative contact rate (RCR)_ of viruses and bacteria within sea ice brines compared to underlying seawater. The relationship between RCR and temperature, calculated by Wells and Deming 2006b, is shown in the repository file [WellsRCR.png](https://github.com/gshowalt/VirusPopModel/blob/main/WellsRCR.png) 

Here's our figure of carbon cycling (1 = maximum carbon) as a function of temperature **WITH** Brine Concentrating Factor

![Practice Text for Sizing](https://github.com/gshowalt/VirusPopModel/blob/main/CE_Temp_RCR_line_withRCRline.jpeg)




### Notes
 <sup>1</sup> While values were parameterized from literature, unmanipulated parameters produced [run-away conditions](https://github.com/gshowalt/VirusPopModel/blob/main/CE_Grid_withRCR_runaway.jpeg) of growth and infection at certain temperatures prior to parameter fitting using the method detail in S1.



### supplemental information

## S1 Parameter fitting

![Fig](https://github.com/gshowalt/VirusPopModel/blob/main/TimeDependent_ParamFit_33_DensityPlot.png)
 
In order to more accurately constrain parameters to the field observations, we wanted to compare calculate VBRs to observed values. As previously mentioned, observations of VBR within sea ice demonstrate high variability, and VBR can reach ratios fo 10,000 : 1 (refer to [Figure 1](https://github.com/gshowalt/VirusPopModel/blob/main/VBRfigure_recreation.png). 

With these steady state solutions, we then wanted to probe for parameter ranges which reproduce observed VBR. To test parameter ranges in VBR, we iterated through random ranges of parameters, calculating VBRs for each set of randomly chosen parameter ranges and comparing the calculated VBR distribution to the observed distributions in the above figure. When the calculated VBR distribution was considered "the same" as the observed VBR distribution (i.e., 95% by a [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test)), we collected that set of parameter values as a possible "solution"  - a set of parameter ranges which could be reflective of the true environment. 46 (I gave up waiting on the program, for the paper we can do more/an even number) of these parameter ranges are plotted below, with the range shown for burst size, growth rate, and decay rate by individual lines and the average value in the circles. The ranges are ordered in ascending order, top to bottom, according to adsorption rate (find the data [here](https://github.com/gshowalt/VirusPopModel/blob/main/TDParams_runs46.xlsx))

The code for this step is given in [This file](https://github.com/gshowalt/VirusPopModel/blob/main/Code/TD_ParamFitting.py)


## S2 Steady state behavior 

The steady state solutions for the system were determined by hand and are given below:

> B = delta/(phi*(beta - 1))

> V = ((mu * N)/(N+Q) - d) / (gamma * phi)

> N = (n * z * d * Q)/(alpha * (g-1) + n * z * (mu - d))


  
