# HePy

[![DOI](https://zenodo.org/badge/593777714.svg)](https://zenodo.org/badge/latestdoi/593777714)

A python code for modelling helium (He) cooling ages in apatite and zircon, following the approach of Ketcham[^1]:
- Designed as a python module, leveraging scipy and other standard libraries
- Includes apatite[^3] and zircon[^5] radiation damage and annealing models (RDAAM)
- RDAAM implementations are based on the DAAM matlab scripts by W. Guenthner [^2]
- Also includes simple diffusion as a function of temperature, in apatite [^4] and zircon [^6]

### Example: Calculating an apatite helium cooling age using the RDAAM diffusion model
```
import numpy as np
from HePy import Model_He

t = np.arange(60,0,-0.5)
T = np.arange(120,0,-4)

ap_1 = Model_He(mineral='apatite',U=10,Th=40,r=100)
print(ap_1.solve(t,T,k='rdaam'))
```
### Convergence of the model age
The model *generally* converges where $dt \leq 0.5  Myrs$, and so tT paths where $\max(dt) \geq 0.5  Myrs$ should be resampled (e.g., using linear interpolation) **before** running the *Model_He.solve()* function. Care should be taken as this is assumed based on testing where $\max (dT/dt < 20 ^{\circ}C/Myr)$ and $T >> Tclosure$. The module will be improved with a resampling function in a future update.

[^1]: Ketcham, R.A., 2005, Forward and inverse modeling of low-temperature thermochronometry data: Reviews in Mineralogy and Geochemistry , v. 58, no. 1, p. 275–314, https://doi.org/10.2138/rmg.2005.58.11
[^2]: Guenthner, W.R., 2020, wrguenthner/DAAM: Second release of damage accumulation and annealing models. (v1.2). Zenodo. https://doi.org/10.5281/zenodo.4289246
[^3]: Flowers, R.M., Ketcham, R.A., Shuster, D.L., and Farley, K.A., 2009, Apatite (U-Th)/He thermochronometry using a radiation damage accumulation and annealing model: Geochimica et Cosmochimica Acta , v. 73, no. 8, p. 2347–2365, https://doi.org/10.1016/j.gca.2009.01.015
[^4]: Farley, K.A., 2000, Helium diffusion from apatite: General behavior as illustrated by Durango fluorapatite: Journal of Geophysical Research. Solid Earth , v. 105, p. 2903–2914, https://doi.org/10.1029/1999JB900348
[^5]: Guenthner, W.R., Reiners, P.W., Ketcham, R.A., Nasdala, L., and Giester, G., 2013, Helium diffusion in natural zircon: Radiation damage, anisotropy, and the interpretation of zircon (U-Th)/He thermochronology: American Journal of Science , v. 313, no. 3, p. 145–198, https://doi.org/10.2475/03.2013.01
[^6]: Reiners, P.W., Spell, T.L., Nicolescu, S., and Zanetti, K.A., 2004, Zircon (U-Th)/He thermochronometry: He diffusion and comparisons with 40Ar/39Ar dating: Geochimica et Cosmochimica Acta, v. 68, no. 8, p. 1857–1887, https://doi.org/10.1016/j.gca.2003.10.021
