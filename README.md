# HePy
A python implementation of 'HeFTy': A software package for modelling He cooling ages in apatite and zircon [^1]

- Based on the DAAM matlab scripts created by Willy Guenther [^4]
- Redesigned as a module and optimised for python
- Leverages scipy and other standard python libraries

### Diffusion models
- Radiation damage and annealing (RDAAM) [^2]
- Zircon???
- Thermally-controlled He diffusion [^3]

### Modelling approach

Add flow diagram here!

### Example: Calculating an apatite helium cooling age using the RDAAM diffusion model
```
import numpy as np
from HePy_v0.1 import Model_He

t = np.arange(60,0,-2)
T = np.arange(120,0,-4)

ap_1 = Model_He(mineral='apatite',U=20,Th=2,r=70)
ap1_age = apatite_1.solve(t,T,k='rdaam')
```
### Numerical stability of the 1D diffusion solution
tT history should be resampled

### How to acknowledge this code
If you use HePy in your scientific work, please acknowledge it as:

> HePy v.01 - A python implementation of HeFTy DOI:

The key references below should also be referenced (depending on the diffusion model used).

[^1]: Ketcham, 2005
[^2]: Flowers et al., 2009
[^3]: Farley et al., 2000
[^4]: DAAM matlab script by Willy Guenthner https://github.com/wrguenthner/DAAM
