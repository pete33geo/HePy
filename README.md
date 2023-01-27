# HePy
A python implementation of the 'HeFTy' modelling approach of Ketcham (2005)

## Theory
HeFTy [^1]
RDAAM [^2]
Thermal Diffusion [^3]
DAAM[^4]

## Modelling approach

HePy leverages scipy to solve the 1D diffusion equation.

### Example: Calculating an apatite helium cooling age using the RDAAM diffusion model
```
import numpy as np
from HePy_v0.1 import Model_He

t = np.arange(60,0,-2)
T = np.arange(120,0,-4)

ap_1 = Model_He(mineral='apatite',U=20,Th=2,r=70)
ap1_age = apatite_1.solve(t,T,k='rdaam')
```
## Numerical stability of the 1D diffusion solution
tT history should be resampled

## How to acknowledge this code
If you use HePy in your scientific work, please acknowledge it

> HePy v.01 - A python implementation of HeFTy DOI:

[^1]: Ketcham, 2005
[^2]: Flowers et al., 2009
[^3]: Farley et al., 2000
[^3]: DAAM matlab script by Willy Guenthner https://github.com/wrguenthner/DAAM
