########################################################################################################################
#
# Mie Scattering program
#
# Nuno de Sousa, MoLE group
# November 2019
#
########################################################################################################################

import numpy as np
np.seterr(divide='ignore', invalid='ignore')

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mole_mie.MieScatt import MieScatt
from mole_mie.Materials import Materials
from mole_mie.MieScatt import MappingSpace

if __name__ == '__main__':

    wavelength_start = 1000  # wavelength range lower limit
    wavelength_end = 2000    # wavelength range upper limit
    wavelength_steps = 500     # number of intervals
    wavelength_list = np.arange(wavelength_start, wavelength_end,
                                (wavelength_end - wavelength_start) / wavelength_steps)

    radius = 230
    medium_ref_index = 1
    N_multipoles = 10

    mie = MieScatt()
    Si = Materials('Si', source='constant', value = 3.5+0j)
    mie.set_params(radius=radius, medium_ref_index=medium_ref_index, N_multipoles=N_multipoles)
    mie.scan_cross_sections(wavelength_list, Si)

    maps = MappingSpace()
    maps.mesh_gen(xmin=-500, xmax=500, xstep=10, ymin=-500, ymax=500, ystep=10, plane='xz')
    maps.convert_cartesian_to_spherical()

    mie.compute_fields(1126, Si, 10, maps.points, components='EH', rat='Scatt')
    mie.abs_fields()

    fig, ax = plt.subplots(figsize=(9, 5.5))

    ax.plot(mie.cross_sections['Qabs'], label='$Q_{abs}$')
    ax.plot(mie.cross_sections['Qext'], label='$Q_{ext}$', linestyle=':')
    ax.plot(mie.cross_sections['Qscat'], label='$Q_{scat}$')
    ax.plot(mie.expanded_Qext['Qext_a1'], label='$Q_{ext, a1}$')
    ax.plot(mie.expanded_Qext['Qext_b1'], label='$Q_{ext, b1}$')
    ax.plot(mie.expanded_Qext['Qext_b2'], label='$Q_{ext, b2}$')
    ax.set_xlabel('$\lambda$ (m)', fontsize=16)
    ax.set_ylabel('', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.tick_params(axis='both', which='minor', labelsize=13)
    plt.legend(fontsize=13)
    plt.title('Cross Sections', fontsize=14)
    plt.xlim(1000E-9, 2000E-9)
    plt.show();

    fig, ax = plt.subplots(figsize=(8, 8))
    z = np.real(
        mie.E_scatt_int_abs.reshape(int(np.sqrt(len(mie.E_scatt_int_abs))), int(np.sqrt(len(mie.E_scatt_int_abs)))))
    im = ax.imshow(z, cmap=cm.jet, interpolation='bilinear', extent=[-500, 500, -500, 500], vmin=0, vmax=9,
                   origin='lower')
    # ax.contour(z, levels = 10, colors='blue', linewidths = 0.5, origin='lower', extent=[-1000,1000,-1000,1000])
    ax.set_xlabel('x', fontsize=16)
    ax.set_ylabel('z', fontsize=16)
    plt.tick_params('both', labelsize=13)
    # fig.colorbar(im, ax = ax)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

    circle = plt.Circle((0, 0), 230, color='black', clip_on=False, fill=False, linewidth=3)
    # circle = plt.Circle((0, 0), 230, color='black', clip_on=False, fill = True)
    ax.add_artist(circle)
    plt.show();

    fig, ax = plt.subplots(figsize=(8, 8))
    z = np.real(
        mie.H_scatt_int_abs.reshape(int(np.sqrt(len(mie.E_scatt_int_abs))), int(np.sqrt(len(mie.E_scatt_int_abs)))))
    im = ax.imshow(z, cmap=cm.jet, interpolation='bilinear', extent=[-500, 500, -500, 500], vmin=0, vmax=9,
                   origin='lower')
    # ax.contour(z, levels = 10, colors='blue', linewidths = 0.5, origin='lower', extent=[-1000,1000,-1000,1000])
    ax.set_xlabel('x', fontsize=16)
    ax.set_ylabel('z', fontsize=16)
    plt.tick_params('both', labelsize=13)
    # fig.colorbar(im, ax = ax)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

    circle = plt.Circle((0, 0), 230, color='black', clip_on=False, fill=False, linewidth=3)
    # circle = plt.Circle((0, 0), 230, color='black', clip_on=False, fill = True)
    ax.add_artist(circle)
    plt.show();