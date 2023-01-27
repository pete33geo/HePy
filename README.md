# HePy
A python implementation of the 'HeFTy' modelling approach of Ketcham (2005)

## Second title

## Flow diagram


### Example: Calculating an apatite helium cooling age using the RDAAM diffusion model
```
import numpy as np
from HePy_v0.1 import Model_He

t = np.arange(60,0,-2)
T = np.arange(120,0,-4)

ap_1 = Model_He(mineral='apatite',U=20,Th=2,r=70)
ap1_age = apatite_1.solve(t,T,k='rdaam')
```
