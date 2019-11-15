# MoLE_Mie: A Mie Scattering calculator

Nuno de Sousa and Juan José Saénz  
Donostia International Physics Center (DIPC), 20018 Donostia-San Sebastian, Spain. 

# How to install

Open a terminal and write:
`pip install git+https://github.com/nunodsousa/MieScattering.git`
The MoLE_Mie and all dependences are installed.

## How to use it - A simple introduction

Let us give an example. We will compute the cross sections of a `Si` particle with `r = 230nm` in the range of `1000` and `2000nm`.
To easy visualize the results we recommend the use of Jupyter notebook.
 
### Step 1: Import libraries 

Start by importing the libraries.

`from mole_mie.MieScatt import MieScatt     # MieScatt library to compute the Mie Solution`   
`from mole_mie.Materials import Materials   # Database with materials`  
`import numpy as np                         # Library to create the wavelength list`  
`import matplotlib.pyplot as plt            # Library to plot the result in the end`  

### Step 2: Define parameters

Definition the parameters for the computation. 

`wavelength_start = 1.0E-6  # wavelength range lower limit`  
`wavelength_end = 2.0E-6    # wavelength range upper limit`  
`wavelength_steps = 500     # Number of points to be plotted`    
`wavelength_list = np.arange(wavelength_start, wavelength_end,
                                (wavelength_end - wavelength_start) / wavelength_steps)`  
`radius = 230E-9            # radius of the particle in meters`   
`medium_ref_index = 1       # refractive index of the hosting medium`  
`N_multipoles = 10          # Number of multipoles to be considered`  

### Step 3: Create the material

Select the material from a pre-loaded list.

`Si = Materials('Si')       # Material of the particle`

### Step 4: Create the object MieScatt and run the computation

Creation of the object, copy the parameters and run a scan in wavelengths. It will return the Absorption, Extintion and Scattering Cross Sections.

`mie = MieScatt()           # Initialization of the object`  
`mie.set_params(radius=radius, medium_ref_index=medium_ref_index, N_multipoles=N_multipoles)`  
`mie.scan_cross_sections(wavelength_list, Si) # scan in wavelengths` 

### Extra Step: Do you want to see the result?

Finally the data is ready for plot. We use the `matplotlib library to execute this step.

`fig, ax = plt.subplots(figsize = (9, 5.5))`  
`ax.plot(mie.cross_sections['Qabs'], label = '$Q_{abs}$')`
`ax.plot(mie.cross_sections['Qext'], label = '$Q_{ext}$', linestyle = ':')`
`ax.plot(mie.cross_sections['Qscat'], label = '$Q_{scat}$')`
`ax.set_xlabel('$\lambda$ (nm)', fontsize = 16)`  
`ax.set_ylabel('', fontsize = 16)`  
`ax.tick_params(axis='both', which='major', labelsize=13)`  
`ax.tick_params(axis='both', which='minor', labelsize=13)`  
`plt.legend(fontsize = 13)`  
`plt.title('Cross Sections', fontsize = 14)`  
`plt.xlim(1000E-9, 2000E-9);`  