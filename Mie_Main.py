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

from mole_mie.MieScatt import MieScatt
from mole_mie.Materials import Materials

if __name__ == '__main__':
    wavelength_start = 1.0E-6  # wavelength range lower limit
    wavelength_end = 2.0E-6    # wavelength range upper limit
    wavelength_steps = 500     # number of intervals
    wavelength_list = np.arange(wavelength_start, wavelength_end,
                                (wavelength_end - wavelength_start) / wavelength_steps)

    radius = 230E-9
    medium_ref_index = 1
    N_multipoles = 10
    Si = Materials('Si')

    mie = MieScatt()
    mie.set_params(radius=radius, medium_ref_index=medium_ref_index, N_multipoles=N_multipoles)
    mie.scan_cross_sections(wavelength_list, Si)

    fig, ax = plt.subplots(figsize=(9, 5.5))

    ax.plot(mie.cross_sections['Qabs'], label='$Q_{abs}$')
    ax.plot(mie.cross_sections['Qext'], label='$Q_{ext}$', linestyle=':')
    ax.plot(mie.cross_sections['Qscat'], label='$Q_{scat}$')
    ax.set_xlabel('$\lambda$ (m)', fontsize=16)
    ax.set_ylabel('', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.tick_params(axis='both', which='minor', labelsize=13)
    plt.legend(fontsize=13)
    plt.title('Cross Sections', fontsize=14)
    plt.xlim(1000E-9, 2000E-9)
    plt.show();