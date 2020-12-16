# VirusPopModel

### Abstract 

Observations indicating high concentrations of viruses within sea-ice brines relative to bacteria indicate that bacteriophage may play an outsized role in shaping host community structure and recycling organic matter within sea ice. However, while these high virus-to-bacteria ratios may be a function of active viral production and release, they may also be the result of slow virion decay due to the low temperature, high salinity conditions of sea-ice brines. Here, we present a mathematical model of virus–host interactions within sea-ice brines, one that uses viral and bacterial abundance data to interpret which mechanisms may most influence population dynamics. We also assess the existence and potential impact of a viral shunt within sea ice. Data from both field samples and laboratory isolates were used to achieve most likely parameter distributions for *in situ* communities, constraining the model to observed dynamics. Deeper understanding of the viral impact on sea-ice communities and the recycling of dissolved organic carbon within the ice will help to understand how these largely heterotrophic communities maintain activity under extreme conditions and further constrain Arctic carbon budgets through the seasons.

# Introduction
This project was written to demonstrate the potential contribution of bacteriophage viruses carbon cycling within sea ice brines. Though within the cold, hypersaline interstitial waters of sea ice microbial activity is reduced compared to more temperate waters, physical concentration may enhance the contact rate between bacteria and their viruses to an extent which allows significant viral infection to take place. 

**We hypothesize that, as a result of physical concentration enhancing viral infection, viruses within sea ice brines may play an out-sized role in facilitating carbon cycling by contributing dissolved organic carbon to sea ice brine DOC pools in the form of lytic material.**

Here, we have built a simple population dynamic model using ordinary differential equations to demonstrate potential and constrain rates of viral infection and virally-mediated carbon cycling within sea ice brines. 


### The system and its challenges
[Observations of VBR within sea ice](https://github.com/gshowalt/VirusPopModel/blob/main/VBRfigure_recreation.png) demonstrate high variability, and VBR can reach ratios fo 10,000 : 1 - much higher than the typical 10:1 or 100:1 values seen in seawater (Figure 1, below).

![Fig](https://github.com/gshowalt/VirusPopModel/blob/main/VBRfigure_recreation.png)


### Introducing the problem
Challenges of understanding




# Methods

### The Equations ###

Equations in this system were modeled on those from given in [*Quantitative Viral Ecology* by Joshua Weitz (2015)](https://press.princeton.edu/books/hardcover/9780691161549/quantitative-viral-ecology) These modles use simple differential equations to build interacting populations of bacteria (B) and viruses (V) by modeling their change as a function of time and nutrient concentration (N), such as:

Bacterial Population
> dBdt = growth - infection - death 

> dBdt = µ * N -  φ * V * B - δb * B 

Virus Population
> dVdt = production (*lysis*) - adsorption - decay 

> dVdt = β * φ * V * B -  1 * φ * V * B - δv * V 


Which were then coded into Python 3:
```
dNdt = (-alpha * (N / (N + Q)) * B) + (g * (alpha  * (N/(N+Q))*B)) + (n * 1e-7 * phi * V * B)
dBdt = ((mu) * (N/(Q + N)) * B) - (phi * V * B) - d*B
dVdt =  (beta * B * phi * V) - (phi * V * B) -  (m * V)
```

Include the RCR as a key element of the model to convey the sea-iciness


### Assumptions ###
These equations differ from those by Weitz or similar open-water models in that they _do not contain any nutrient inflow_ - instead, they consider the sea ice brine pore as a closed system. As a result - in the equations given above - nutrients are consumed in the given equations and only recycled through exudate and lysis (not through bacterial death or viral decay).
### Parameters ###
Loosey bound parameters from literature - give the spreadsheet here

### "Experiments ###


# Results & Discussion
## 1. Can we replicate Virus to Bactera Ratios by tuning parameters?

Collecting parameter values from literature produced a set of ranges given in [Table 1](https://github.com/gshowalt/VirusPopModel/blob/main/Table1.png) (do we want to include this?).


When we run the experiments, we see either extinction of populations or runaway growth of virus populations, as shown below in Figure S1.

![Fig](https://github.com/gshowalt/VirusPopModel/blob/main/TimeDependent_withRCRnoManipulation.png)
 
This indicates that the parameters given in literature may not be accurate for the given system (or the equations themselves are missing a vital element of the system).



In order to more accurately constrain parameters to the field observations, we wanted to compare calculate VBRs to observed values. As previously mentioned, observations of VBR within sea ice demonstrate high variability, and VBR can reach ratios fo 10,000 : 1 (refer to [Figure 1](https://github.com/gshowalt/VirusPopModel/blob/main/VBRfigure_recreation.png). 


The steady state solutions for the system were determined by hand and are given below:

> B = detla/(phi*(beta - 1))

> V = ((mu * N)/(N+Q) - d) / (gamma * phi)

> N = (n * z * d * Q)/(alpha * (g-1) + n * z * (mu - d))

With these steady state solutions, we then wanted to probe for parameter ranges which reproduce observed VBR. To test parameter ranges in VBR, we iterated through random ranges of parameters, calculating VBRs for each set of randomly chosen parameter ranges and comparing the calculated VBR distribution to the observed distributions in the above figure. When the calculated VBR distribution was considered "the same" as the observed VBR distribution (i.e., 95% by a [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test)), we collected that set of parameter values as a possible "solution"  - a set of parameter ranges which could be reflective of the true environment. 46 (I gave up waiting on the program, for the paper we can do more/an even number) of these parameter ranges are plotted below, with the range shown for burst size, growth rate, and decay rate by individual lines and the average value in the circles. The ranges are ordered in ascending order, top to bottom, according to adsorption rate (find the data [here](https://github.com/gshowalt/VirusPopModel/blob/main/TDParams_runs46.xlsx))

![Parameter Distributions](https://github.com/gshowalt/VirusPopModel/blob/main/TimeDependent_ParamFit_46.png)
The code for this step is given in [This file](https://github.com/gshowalt/VirusPopModel/blob/main/TD_ParamFitting.py)

**We are also re-running this with absolute abundance fitting to narrow this range, which will hopefully allow us to better reproduce time-dependent dynamics**

## 2. Do we want to explore lysis vs. lysogeny?


## 3. How does viral infection contribute to carbon cycling within sea ice?
First, we use values collected from literature to show that without physical concentration due to _brine concentrating factor_ (BCF), viral infection would have negligble impact on microbial populations and carbon cycling within sea ice.

In order to make this demonstration, we first took our above equations coded into Python and applied biological values (i.e. beta (burst size), phi (adsorption rate), mu (bacterial growth rate), and delta (viral decay rate)) collected from literature and parameterized as function dependent on temperature  <sup>[1](###Notes)</sup>. Temperature-dependent values allowed us to 

The code for this step is uploaded to the repository as [CarbonEquiv_noRCR.py](https://github.com/gshowalt/VirusPopModel/edit/main/CarbonEquiv_Talmy.py) .

Here's our figure of carbon cycling (1 = maximum carbon) as a function of temperature **WITHOUT** Brine Concentrating Factor

![Practice Text for Sizing](https://github.com/gshowalt/VirusPopModel/blob/main/CE_Temp_noRCR_line.jpeg)






Secondly, we apply a **physical concentration parameter** to the model. This physical concentration, a result of tightly constricted pore space within sea ice crystals, has been suggested to increase the _relative contact rate (RCR)_ of viruses and bacteria within sea ice brines compared to underlying seawater. The relationship between RCR and temperature, calculated by Wells and Deming 2006b, is shown in the repository file [WellsRCR.png](https://github.com/gshowalt/VirusPopModel/blob/main/WellsRCR.png) 

Here's our figure of carbon cycling (1 = maximum carbon) as a function of temperature **WITH** Brine Concentrating Factor

![Practice Text for Sizing](https://github.com/gshowalt/VirusPopModel/blob/main/CE_Temp_RCR_line_withRCRline.jpeg)



# Additional Considerations

Adsorption
Salt precipitation


### Notes
 <sup>1</sup> While values were parameterized from literature, unmanipulated parameters produced [run-away conditions](https://github.com/gshowalt/VirusPopModel/blob/main/CE_Grid_withRCR_runaway.jpeg) of growth and infection at certain temperatures. 



#### notes for self (Max)
Follow up from 11/30 mtg --> postponed to Thursday
DONE 1. create time-dependent plot - > does it behave as expected? Do we see a net loss of carbon due to *mV* term
2. Check on carbon cycling equations - do they still work if you pull them out of the odeint term
DONE 3. make carbon cycling plot w/o breakdown into lysate/virions
4. think of most intuitive figure to show RCR relationship in carbon cycling figure
DONE 5. compare parameter fitting for simple vs. complex eqs.
DONE 6. craft narrative order of the story
