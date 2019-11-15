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
import pandas as pd

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

    __version__ = '0.1.1'
    __codename__ = "ingrown toenail"

    a = None
    m = None
    mt = None
    N_multipoles = None

    coeffs = None
    expanded_Qext = None
    cross_sections = None

    def __init__(self):
        self.showinfo()

    def showinfo(self):
        """

        :return:
        """
        print(colored(14 * '_______', 'blue'))
        print(colored('Mie Scattering program Initiated.\n', 'blue'))
        print('You are using {}. Please cite the following:'.format(colored('MoLE_Mie', 'red')))
        print('* N. de Sousa, J.J. Saénz, \"The title of the paper\".\n')
        print('If you use the pre-loaded database, please cite:')
        print('* M.N.Polyanskiy, \"Refractive index database,\" https: // refractiveindex.info. Accessed on 2019-10-04.')
        print('Version', self.__version__, self.__codename__)
        print(colored(14*'_______', 'blue'))

    def set_params(self, radius=None, medium_ref_index=None, N_multipoles=None):
        """

        :param radius:
        :param medium_ref_index:
        :param N_multipoles:
        :return:
        """
        self.a = radius;
        self.m = medium_ref_index;
        self.N_multipoles = N_multipoles;
        print("Loaded Parameters for Si")
        print("radius = ", self.a, "(m)")
        print("refractive index of the medium = ", self.m)
        print("Number of multipoles = ", self.N_multipoles)

    def compute_coeffs(self, mt, m, mp, alpha, beta):
        """

        :param mt:
        :param m:
        :param mp:
        :param alpha:
        :param beta:
        :return:
        """
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
        """

        :param wavelength_list:
        :param material:
        :return:
        """

        self.check_parameters()

        N_multipoles = self.N_multipoles;
        m = self.m;
        a = self.a;

        mp = np.arange(1, N_multipoles + 1, 1)  # multipole list

        result = []
        result_Qexta_expanded = []
        result_Qextb_expanded = []

        result_coeff_an = []
        result_coeff_bn = []
        result_coeff_cn = []
        result_coeff_dn = []


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
            
            Qext_a = (2 / alpha ** 2) * (2 * mp + 1) * (an).real
            Qext_b = (2 / alpha ** 2) * (2 * mp + 1) * (bn).real

            result.append((wavelength, Qscat, Qext))
            result_Qexta_expanded.append(Qext_a)
            result_Qextb_expanded.append(Qext_b)

            result_coeff_an.append(an)
            result_coeff_bn.append(bn)
            result_coeff_cn.append(cn)
            result_coeff_dn.append(dn)

        result = np.array(result)

        col_names_a = ['Qext_a' + str(i) for i in range(1,N_multipoles + 1)]
        col_names_b = ['Qext_b' + str(i) for i in range(1, N_multipoles + 1)]

        col_coeff_a = ['a' + str(i) for i in range(1, N_multipoles + 1)]
        col_coeff_b = ['b' + str(i) for i in range(1, N_multipoles + 1)]
        col_coeff_c = ['c' + str(i) for i in range(1, N_multipoles + 1)]
        col_coeff_d = ['d' + str(i) for i in range(1, N_multipoles + 1)]

        result_Qexta_expanded = pd.DataFrame(np.array(result_Qexta_expanded), columns=col_names_a)
        result_Qextb_expanded = pd.DataFrame(np.array(result_Qextb_expanded), columns=col_names_b)

        result_coeff_an = pd.DataFrame(np.array(result_coeff_an), columns=col_coeff_a)
        result_coeff_bn = pd.DataFrame(np.array(result_coeff_bn), columns=col_coeff_b)
        result_coeff_cn = pd.DataFrame(np.array(result_coeff_cn), columns=col_coeff_c)
        result_coeff_dn = pd.DataFrame(np.array(result_coeff_dn), columns=col_coeff_d)

        Qscat = np.sum(np.stack(list(result[:, 1])), axis=1)
        Qext = np.sum(np.stack(list(result[:, 2])), axis=1)
        Qabs = Qext - Qscat

        self.Qscat = Qscat
        self.Qext = Qext
        self.Qabs = Qabs

        result_Qexta_expanded.index = wavelength_list
        result_Qextb_expanded.index = wavelength_list

        self.expanded_Qext = pd.concat([result_Qexta_expanded, result_Qextb_expanded], axis = 1, sort = False)

        self.coeffs = pd.concat([result_coeff_an, result_coeff_bn, result_coeff_cn, result_coeff_dn], axis = 1, sort = False)
        self.coeffs.index = wavelength_list

        self.cross_sections = pd.DataFrame([self.Qabs, self.Qext, self.Qscat]).T.rename(columns = {0:'Qabs', 1:'Qext', 2: 'Qscat'})
        self.cross_sections.index = wavelength_list

    def check_parameters(self, *arg):
        """

        :param arg:
        :return:
        """
        if (self.a == None):
            raise ValueError('The particle Radius is not defined.')
        if (self.m == None):
            raise ValueError('The medium refractive index is not defined.')
        if (self.N_multipoles == None):
            raise ValueError('Number of multipoles to be considered is not defined.')

        for par in arg:
            if (par == None):
                raise ValueError