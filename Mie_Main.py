########################################################################################################################
#
# Mie Scattering program
#
# Nuno de Sousa, MoLE group
# November 2019
#
########################################################################################################################

import numpy as np
import matplotlib.pyplot as plt

from libs.MieLib import MieScatt
from libs.Materials import Materials

if __name__ == '__main__':
    wavelength_start = 1.0E-6  # wavelength range lower limit
    wavelength_end = 2.0E-6  # wavelength range upper limit
    wavelength_steps = 500
    wavelength_list = np.arange(wavelength_start, wavelength_end,
                                (wavelength_end - wavelength_start) / wavelength_steps)

    radius = 230E-9
    medium_ref_index = 1
    N_multipoles = 10
    Si = Materials('Si')

    mie = MieScatt()
    mie.set_params(radius=radius, medium_ref_index=medium_ref_index, N_multipoles=N_multipoles)
    mie.scan_cross_sections(wavelength_list, Si)

    # other way of call the cross section
    #mie.cross_sections(1260E-9, Si)

    fig, ax = plt.subplots(figsize = (9, 5.5))

    ax.plot(wavelength_list*1E9, mie.Qext, label = '$Q_{ext}$')
    ax.plot(wavelength_list*1E9, mie.Qscat, label = '$Q_{scat}$', linestyle = ':')
    ax.plot(wavelength_list*1E9, mie.Qabs, label = '$Q_{abs}$')
    ax.set_xscale('log')
    ax.set_xlabel('$\lambda$ (nm)', fontsize = 16)
    ax.set_ylabel('', fontsize = 16)
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.tick_params(axis='both', which='minor', labelsize=13)
    plt.legend(fontsize = 13)
    plt.title('Cross Sections', fontsize = 14)
    plt.xlim(1000, 2000);
    plt.show();