# MieScattering
Mie Scattering calculator.  
Nuno de Sousa and Juan José Saénz  
Donostia International Physics Center.  

# How to install

Open a terminal and write:
`pip install git+https://github.com/nunodsousa/MieScattering.git`
And it should be ready.

## How to use it

The simplest way is to give an example. We will compute the cross sections of a Si particle with r = 230nm in the range of E-6 and 2E-6m.
 
To start, probably the best solution is to use jupyter notebook, so you can start by open it.
### Step 1: Import libraries 

`from mole_mie.MieScatt import MieScatt`  
`from mole_mie.Materials import Materials`

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

### Extra Step: Do you want to see the result?
`fig, ax = plt.subplots(figsize = (9, 5.5))`  
`ax.plot(wavelength_list*1E9, mie.Qext, label = '$Q_{ext}$')`  
`ax.plot(wavelength_list*1E9, mie.Qscat, label = '$Q_{scat}$', linestyle = ':')`  
`ax.plot(wavelength_list*1E9, mie.Qabs, label = '$Q_{abs}$')`  
`ax.set_xscale('log')`  
`ax.set_xlabel('$\lambda$ (nm)', fontsize = 16)`  
`ax.set_ylabel('', fontsize = 16)`  
`ax.tick_params(axis='both', which='major', labelsize=13)`  
`ax.tick_params(axis='both', which='minor', labelsize=13)`  
`plt.legend(fontsize = 13)`  
`plt.title('Cross Sections', fontsize = 14)`  
`plt.xlim(1000, 2000);`  
`plt.show();`