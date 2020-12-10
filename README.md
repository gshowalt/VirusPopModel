# VirusPopModel

### quick abstract here
This project was written to demonstrate the potential contribution of bacteriophage viruses carbon cycling within sea ice brines. Though within the cold, hypersaline interstitial waters of sea ice microbial activity is reduced compared to more temperate waters, physical concentration may enhance the contact rate between bacteria and their viruses to an extent which allows significant viral infection to take place. 

**We hypothesize that, as a result of physical concentration enhancing viral infection, viruses within sea ice brines may play an out-sized role in facilitating carbon cycling by contributing dissolved organic carbon to sea ice brine DOC pools in the form of lytic material.**

Here, we have built a simple population dynamic model using ordinary differential equations to demonstrate potential and constrain rates of viral infection and virally-mediated carbon cycling within sea ice brines. 

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

These equations differ from those by Weitz or similar open-water models in that they _do not contain any nutrient inflow_ - instead, they consider the sea ice brine pore as a closed system. As a result - in the equations given above - nutrients are consumed in the given equations and only recycled through exudate and lysis (not through bacterial death or viral decay).


## 1. Can we replicate Virus to Bactera Ratios by tuning parameters?

In addition to understanding virally mediated carbon flow within sea ice brines, we want to undertand potential controls on the virus to bacteria (VBR) ratio in sea ice. [Observations of VBR within sea ice](https://github.com/gshowalt/VirusPopModel/blob/main/VBRfigure_recreation.png) demonstrate high variability, and VBR can reach ratios fo 10,000 : 1 - much higher than the typical 10:1 or 100:1 values seen in seawater.

![Fig](https://github.com/gshowalt/VirusPopModel/blob/main/VBRfigure_recreation.png)

To test parameter ranges in VBR, we iterated through random ranges of parameters in the simplest possible system (insert EQs), calculating VBRs for each set of randomly chosen parameter ranges and comparing the calculated VBR distribution to the observed distributions in the above figure. When the calculated VBR distribution was considered "the same" as the observed VBR distribution (i.e., 95% by a [Kolmogorov-Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov–Smirnov_test)), we collected that set of parameter values as a possible "solution"  - a set of parameter ranges which could be reflective of the true environment. 50 of these parameter ranges are plotted below, with the range shown for burst size, growth rate, and decay rate by individual lines and the average value in the circles. The ranges are ordered in ascending order, top to bottom, according to decay rate average.

![Parameter Distributions](https://github.com/gshowalt/VirusPopModel/blob/main/ParameterDistribution.png)

From this figure, we can see that growth rate and decay rate appear to tightly co-vary, while burst size changs independent to either growth rate or burst size. Additionally, the wide range of burst size values for each parameter set implies that burst size has minimal effect on the final VBR.

The code for this step is given in [This file](https://github.com/gshowalt/VirusPopModel/blob/main/ParameterFitting.py)

## 2. Can we replicate absoute abundance of viruses and bacteria by tuning parameters?

We repeated the process from above against absolute abundance values, finding parameters that can replicate absolute abundance for both viruses and bacteria within KS metric > 0.05. 

![FittoBoth](https://github.com/gshowalt/VirusPopModel/blob/main/FittoBoth.png)


## 3. How does viral infection contribute to carbon cycling within sea ice?
First, we use values collected from literature to show that without physical concentration due to _brine concentrating factor_ (BCF), viral infection would have negligble impact on microbial populations and carbon cycling within sea ice.

In order to make this demonstration, we first took our above equations coded into Python and applied biological values (i.e. beta (burst size), phi (adsorption rate), mu (bacterial growth rate), and delta (viral decay rate)) collected from literature and parameterized as function dependent on temperature  <sup>[1](###Notes)</sup>. Temperature-dependent values allowed us to 

The code for this step is uploaded to the repository as [CarbonEquiv_noRCR.py](https://github.com/gshowalt/VirusPopModel/edit/main/CarbonEquiv_Talmy.py) .

Here's our figure of carbon cycling (1 = maximum carbon) as a function of temperature **WITHOUT** Brine Concentrating Factor

![Practice Text for Sizing](https://github.com/gshowalt/VirusPopModel/blob/main/CE_Temp_noRCR_line.jpeg)






Secondly, we apply a **physical concentration parameter** to the model. This physical concentration, a result of tightly constricted pore space within sea ice crystals, has been suggested to increase the _relative contact rate (RCR)_ of viruses and bacteria within sea ice brines compared to underlying seawater. The relationship between RCR and temperature, calculated by Wells and Deming 2006b, is shown in the repository file [WellsRCR.png](https://github.com/gshowalt/VirusPopModel/blob/main/WellsRCR.png) 

Here's our figure of carbon cycling (1 = maximum carbon) as a function of temperature **WITH** Brine Concentrating Factor

![Practice Text for Sizing](https://github.com/gshowalt/VirusPopModel/blob/main/CE_Temp_withRCR_line.jpeg)

### Notes
 <sup>1</sup> While values were parameterized from literature, unmanipulated parameters produced [run-away conditions](https://github.com/gshowalt/VirusPopModel/blob/main/CE_Grid_withRCR_runaway.jpeg) of growth and infection at certain temperatures. 



#### notes for self (Max)
Follow up from 11/30 mtg:
1. create time-dependent plot - > does it behave as expected? Do we see a net loss of carbon due to *mV* term
2. Check on carbon cycling equations - do they still work if you pull them out of the odeint term
3. make carbon cycling plot w/o breakdown into lysate/virions
4. think of most intuitive figure to show RCR relationship in carbon cycling figure
5. compare parameter fitting for simple vs. complex eqs.
6. craft narrative order of the story
