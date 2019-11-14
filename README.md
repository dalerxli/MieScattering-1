# MieScattering
A simple Mie scattering program to compute the exact solution.


## How to use it

### Step 1: Import libraries 

from libs.MieLib import MieScatt
from libs.Materials import Materials

### Step 2: Generate parameters

`wavelength_start = 1.0E-6  # wavelength range lower limit`  
`wavelength_end = 2.0E-6  # wavelength range upper limit`  
`wavelength_steps = 500`  
`wavelength_list = np.arange(wavelength_start, wavelength_end,
                                (wavelength_end - wavelength_start) / wavelength_steps)`  
`radius = 230E-9`   
`medium_ref_index = 1`  
`N_multipoles = 10`  

### Step 3: Create the material

`Si = Materials('Si')`

### Step 4: Create the object MieScatt and run the computation

`mie = MieScatt()`  
`mie.set_params(radius=radius, medium_ref_index=medium_ref_index, N_multipoles=N_multipoles)`  
`mie.scan_cross_sections(wavelength_list, Si)` 