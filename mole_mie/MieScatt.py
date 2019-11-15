###########################################################################################
#
# Mie Library
#
# Nuno de Sousa
# December 2019
###########################################################################################

import numpy as np
import scipy.special as sps
from termcolor import colored

class SpecialFunctions(object):

    def psi(self, n, x):
        return x*sps.spherical_jn(n, x, 0)

    def diff_psi(self, n, x):
        return sps.spherical_jn(n, x, 0) + x*sps.spherical_jn(n, x, 1)

    def xi(self, n, x):
        return x*(sps.spherical_jn(n, x, 0) + 1j*sps.spherical_yn(n, x, 0))

    def diff_xi(self, n,x):
        return (sps.spherical_jn(n, x, 0) + 1j*sps.spherical_yn(n, x, 0)) + x*(sps.spherical_jn(n, x, 1) + 1j*sps.spherical_yn(n, x, 1))

class MieScatt(SpecialFunctions):

    __version__ = '0.1'
    __codename__ = "ingrown toenail"

    a = None
    m = None
    mt = None
    N_multipoles = None

    def __init__(self):
        self.showinfo()

    def showinfo(self):
        print(colored(50 * '__', 'blue'))
        print(colored('Mie Scattering program Initiated.\n', 'blue'))
        print('You are using {}. Please cite the following:'.format(colored('Mie-MoLE', 'red')))
        print('* N. de Sousa, J.J. Saénz, \"The title of the paper\".\n')
        print('If you use the pre-loaded database, please cite:')
        print('* M.N.Polyanskiy, \"Refractive index database,\" https: // refractiveindex.info. Accessed on 2019-10-04.')
        print('Version', self.__version__, self.__codename__)
        print(colored(50*'__', 'blue'))

    def load_Si_particle(self):
        self.a = 0.23E-6;
        self.m = 1;
        self.mt = 3.5 + 0.0j;
        self.N_multipoles = 10;
        print("Loaded Parameters for Si")
        print("radius = ", self.a, "(m)")
        print("refractive index of the medium = ", self.m)
        print("refractive index of the scatterer = ", self.mt)
        print("Number of multipoles = ", self.N_multipoles)

    def set_params(self, radius=None, medium_ref_index=None, N_multipoles=None):
        self.a = radius;
        self.m = medium_ref_index;
        self.N_multipoles = N_multipoles;
        print("Loaded Parameters for Si")
        print("radius = ", self.a, "(m)")
        print("refractive index of the medium = ", self.m)
        print("Number of multipoles = ", self.N_multipoles)

    def compute_coeffs(self, mt, m, mp, alpha, beta):
        an = (mt * self.diff_psi(mp, alpha) * self.psi(mp, beta) - m * self.psi(mp, alpha) * self.diff_psi(mp,
                                                                                                           beta)) / (
                         mt * self.diff_xi(mp, alpha) * self.psi(mp, beta) - m * self.xi(mp, alpha) * self.diff_psi(
                     mp, beta))
        bn = (mt * self.psi(mp, alpha) * self.diff_psi(mp, beta) - m * self.diff_psi(mp, alpha) * self.psi(mp,
                                                                                                           beta)) / (
                         mt * self.xi(mp, alpha) * self.diff_psi(mp, beta) - m * self.diff_xi(mp, alpha) * self.psi(
                     mp, beta))
        cn = (mt * self.xi(mp, alpha) * self.diff_psi(mp, alpha) - mt * self.diff_xi(mp, alpha) * self.psi(mp,
                                                                                                           alpha)) / (
                         mt * self.xi(mp, alpha) * self.diff_psi(mp, beta) - m * self.diff_xi(mp, alpha) * self.psi(
                     mp, beta))
        dn = (mt * self.diff_xi(mp, alpha) * self.psi(mp, alpha) - mt * self.xi(mp, alpha) * self.diff_psi(mp,
                                                                                                           alpha)) / (
                         mt * self.diff_xi(mp, alpha) * self.psi(mp, beta) - m * self.xi(mp, alpha) * self.diff_psi(
                     mp, beta))
        return an, bn, cn, dn

    def scan_cross_sections(self, wavelength_list, material):

        self.check_parameters()

        N_multipoles = self.N_multipoles;
        m = self.m;
        a = self.a;

        mp = np.arange(1, N_multipoles + 1, 1)  # multipole list

        result = []

        for wavelength in wavelength_list:
            mt = material.refractive_index(wavelength);

            k = 2 * np.pi * m / wavelength
            kt = 2 * np.pi * mt / wavelength

            alpha = k * a
            beta = kt * a

            an, bn, cn, dn = self.compute_coeffs(mt, m, mp, alpha, beta)

            # Evaluation of the cross sections
            Qscat = (2 / alpha ** 2) * (2 * mp + 1) * (np.abs(an) ** 2 + np.abs(bn) ** 2)
            Qext = (2 / alpha ** 2) * (2 * mp + 1) * (an + bn).real

            result.append((wavelength, Qscat, Qext))

        result = np.array(result)

        Qscat = np.sum(np.stack(list(result[:, 1])), axis=1)
        Qext = np.sum(np.stack(list(result[:, 2])), axis=1)
        Qabs = Qext - Qscat

        self.Qscat = Qscat
        self.Qext = Qext
        self.Qabs = Qabs

    def cross_sections(self, wavelength, material):
        """
        Return the cross sections decomposed into the fundamental elements of the multipoles.

        """

        self.check_parameters()

        N_multipoles = self.N_multipoles;
        m = self.m;
        a = self.a;
        mp = np.arange(1, N_multipoles + 1, 1)  # multipole list
        mt = material.refractive_index(wavelength);

        k = 2 * np.pi * m / wavelength
        kt = 2 * np.pi * mt / wavelength

        alpha = k * a
        beta = kt * a

        an, bn, cn, dn = self.compute_coeffs(mt, m, mp, alpha, beta)
        # Evaluation of the cross sections
        Qscat = (2 / alpha ** 2) * (2 * mp + 1) * (np.abs(an) ** 2 + np.abs(bn) ** 2)
        Qext = (2 / alpha ** 2) * (2 * mp + 1) * (an + bn).real

        Qscat = np.stack(Qscat)
        Qext = np.stack(Qext)
        Qabs = Qext - Qscat

        return Qabs, Qext, Qscat

    def check_parameters(self, *arg):
        if (self.a == None):
            raise ValueError('The particle Radius is not defined.')
        if (self.m == None):
            raise ValueError('The medium refractive index is not defined.')
        if (self.N_multipoles == None):
            raise ValueError('Number of multipoles to be considered is not defined.')

        for par in arg:
            if (par == None):
                raise ValueError