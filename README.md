# VirusPopModel

This project was written to demonstrate the potential contribution of bacteriophage viruses carbon cycling within sea ice brines. Though within the cold, hypersaline interstitial waters of sea ice microbial activity is reduced compared to more temperate waters, physical concentration may enhance the contact rate between bacteria and their viruses to an extent which allows significant viral infection to take place. 

**We hypothesize that, as a result of physical concentration enhancing viral infection, viruses within sea ice brines may play an out-sized role in facilitating carbon cycling by contributing dissolved organic carbon to sea ice brine DOC pools in the form of lytic material.**

Here, we have built a simple population dynamic model using ordinary differential equations to demonstrate potential and constrain rates of viral infection and virally-mediated carbon cycling within sea ice brines. 

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
dBdt = ((mu) * (N/(Q + N)) * B) - (phi * V * B - d)
dVdt =  (beta * B * phi * V) - (phi * V * B) -  (m * V)
```

### Running without manipulation 
First, we use values collected from literature to show that without physical concentration due to _brine concentrating factor_ (BCF), viral infection would have negligble impact on microbial populations and carbon cycling within sea ice.

In order to make this demonstration, we first took our above equations coded into Python and applied biological values (i.e. beta (burst size), phi (adsorption rate), mu (bacterial growth rate), and delta (viral decay rate)) collected from literature and parameterized as function dependent on temperature  <sup>[1](###Notes)</sup>. Temperature-dependent values allowed us to 

The code for this step is uploaded to the repository as [CarbonEquiv_noRCR.py](https://github.com/gshowalt/VirusPopModel/blob/main/CarbonEquiv_noRCR.py) .

Here's our figure of carbon cycling (1 = maximum carbon) as a function of temperature **WITHOUT** Brine Concentrating Factor

![Practice Text for Sizing](https://github.com/gshowalt/VirusPopModel/blob/main/CE_Temp_noRCR_line.jpeg)






Secondly, we apply a **physical concentration parameter** to the model. This physical concentration, a result of tightly constricted pore space within sea ice crystals, has been suggested to increase the _relative contact rate (RCR)_ of viruses and bacteria within sea ice brines compared to underlying seawater. The relationship between RCR and temperature, calculated by Wells and Deming 2006b, is shown in the repository file [WellsRCR.png](https://github.com/gshowalt/VirusPopModel/blob/main/WellsRCR.png) 

Here's our figure of carbon cycling (1 = maximum carbon) as a function of temperature **WITH** Brine Concentrating Factor

![Practice Text for Sizing](https://github.com/gshowalt/VirusPopModel/blob/main/CE_Temp_withRCR_line.jpeg)


### Notes

