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
    """
    This class implements the vector spherical Harmonics, as well as the Pi_n and Tau_n function.
    We must Merge this class with the one described in MieScatt.py

    Look to this code:
    http://pymiescatt.readthedocs.io/en/latest/forward.html#MiePiTau

    Nuno de Sousa
    """
    __version__ = 0.1

    def psi(self, n, x):
        return x * sps.spherical_jn(n, x, 0)

    def diff_psi(self, n, x):
        return sps.spherical_jn(n, x, 0) + x * sps.spherical_jn(n, x, 1)

    def xi(self, n, x):
        return x * (sps.spherical_jn(n, x, 0) + 1j * sps.spherical_yn(n, x, 0))

    def diff_xi(self, n, x):
        return (sps.spherical_jn(n, x, 0) + 1j * sps.spherical_yn(n, x, 0)) + x * (
                    sps.spherical_jn(n, x, 1) + 1j * sps.spherical_yn(n, x, 1))

    def Pi_n(self, n, theta):
        mu = np.cos(theta)
        p = np.zeros(int(n))
        if (n > 1):
            p[0] = 1
            p[1] = 3 * mu
            for j in range(2, int(n)):
                p[j] = ((2 * j + 1) * (mu * p[j - 1]) - (j + 1) * p[j - 2]) / j
            return p[-1]
        elif (n == 1):
            return 1
        else:
            raise ValueError('You are defining a n smaller than 1 in the Pi_n function. Pi_n is only\
                                 defined for integer numbers >= 1.')

    def Tau_n(self, n, theta):
        mu = np.cos(theta)
        if (n > 1):
            p = np.zeros(int(n))
            t = np.zeros(int(n))
            p[0] = 1
            p[1] = 3 * mu
            t[0] = mu
            t[1] = 3.0 * np.cos(2 * theta)
            for n in range(2, int(n)):
                p[n] = ((2 * n + 1) * (mu * p[n - 1]) - (n + 1) * p[n - 2]) / n
                t[n] = (n + 1) * mu * p[n] - (n + 2) * p[n - 1]
            return t[-1]
        elif (n == 1):
            return mu
        else:
            raise ValueError('You are defining a n smaller than 1 in the Tau_n function. Tau_n is only\
                             defined for integer numbers >= 1.')

    def Mo1n(self, n, sphere_vector, k, inside=False):
        r, theta, phi = sphere_vector[0], sphere_vector[1], sphere_vector[2]
        return (0,
                np.cos(phi) * self.Pi_n(n, theta) * self.radial_function(n, k, r, inside),
                -np.sin(phi) * self.Tau_n(n, theta) * self.radial_function(n, k, r, inside))

    def Me1n(self, n, sphere_vector, k, inside=False):
        r, theta, phi = sphere_vector[0], sphere_vector[1], sphere_vector[2]
        return (0,
                -np.sin(phi) * self.Pi_n(n, theta) * self.radial_function(n, k, r, inside),
                -np.cos(phi) * self.Tau_n(n, theta) * self.radial_function(n, k, r, inside))

    def No1n(self, n, sphere_vector, k, inside=False):
        r, theta, phi = sphere_vector[0], sphere_vector[1], sphere_vector[2]
        #print("a = ", self.radial_function(n, k, r, inside) / (k * r))
        #print("b = ", self.Diff_rho_sf_rho(n, k, r, inside) / (k * r))
        return (np.sin(phi) * n * (n + 1) * np.sin(theta) * self.Pi_n(n, theta) * self.radial_function(n, k, r, inside) / (k * r),
                np.sin(phi) * self.Tau_n(n, theta) * self.Diff_rho_sf_rho(n, k, r, inside) / (k * r),
                np.cos(phi) * self.Pi_n(n, theta) * self.Diff_rho_sf_rho(n, k, r, inside) / (k * r))

    def Ne1n(self, n, sphere_vector, k, inside=False):
        r, theta, phi = sphere_vector[0], sphere_vector[1], sphere_vector[2]
        return (
        np.cos(phi) * n * (n + 1) * np.sin(theta) * self.Pi_n(n, theta) * self.radial_function(n, k, r, inside) / (k * r),
        np.cos(phi) * self.Tau_n(n, theta) * self.Diff_rho_sf_rho(n, k, r, inside) / (k * r),
        -np.sin(phi) * self.Pi_n(n, theta) * self.Diff_rho_sf_rho(n, k, r, inside) / (k * r))

    def radial_function(self, n, k, r, inside=False):
        if (inside == True):
            return sps.spherical_jn(n, k * r, 0)
        else:
            return self.spherical_hn1(n, k * r)

    def spherical_hn1(self, n, z, derivative=False):
        """ Spherical Hankel Function of the First Kind """
        return sps.spherical_jn(n, z, derivative=False) + 1j * sps.spherical_yn(n, z, derivative=False)

    def Diff_rho_sf_rho(self, n, k, r, inside=False):
        if (inside == True):
            return sps.spherical_jn(n, k * r, 0) + k * r * (-sps.spherical_jn(n, k * r) / (2 * k * r) + 1 / 2 * (sps.spherical_jn(n - 1, k * r) - sps.spherical_jn(n + 1, k * r)))
        else:
            return self.spherical_hn1(n, k * r, 0) + k * r * (-self.spherical_hn1(n, k * r) / (2 * k * r) + 1 / 2 * (self.spherical_hn1(n - 1, k * r) - self.spherical_hn1(n + 1, k * r)))


class MieScatt(SpecialFunctions):
    __version__ = '0.1.4.8'
    __codename__ = "dark ingrown toenail"

    a = None
    m = None
    mt = None
    N_multipoles = None

    coeffs = None
    expanded_Qext = None
    cross_sections = None

    c_const = 1; #(nm/s)

    def __init__(self, silence=False):
        """
        Possiblity of hide the initial message.
        :param silence: True or False
        """
        if (silence == False):
            self.showinfo()

    def showinfo(self):
        """
        :return: Information about the library
        """
        print(colored(14 * '_______', 'blue'))
        print(colored('Mie Scattering program Initiated.\n', 'blue'))
        print('You are using {}. Please cite the following:'.format(colored('MoLE_Mie', 'red')))
        print('* N. de Sousa, J.J. Saénz, \"The title of the paper\".\n')
        print('If you use the pre-loaded database, please cite:')
        print(
            '* M.N.Polyanskiy, \"Refractive index database,\" https: // refractiveindex.info. Accessed on 2019-10-04.')
        print('Version', self.__version__, self.__codename__)
        print(colored(14 * '_______', 'blue'))

    def set_params(self, radius=None, medium_ref_index=None, N_multipoles=None):
        """
        Definition of the parameters. It is self-explained
        :param radius:
        :param medium_ref_index:
        :param N_multipoles:
        :return: None
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
        :param mt: refractive index of the particle
        :param m: refractive index of the medium
        :param mp: multipoles vector
        :param alpha: wavenumber in the medium multiplied by the particle radius
        :param beta: wavenumber in the particle multiplied by the particle radius
        :return: mie coefficients
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
        Computes the Mie Solution for a wavelength list.
        :param wavelength_list:
        :param material:
        :return: Variables in the class. cross_sections, expanded_Qext, coeffs
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

        col_names_a = ['Qext_a' + str(i) for i in range(1, N_multipoles + 1)]
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

        self.expanded_Qext = pd.concat([result_Qexta_expanded, result_Qextb_expanded], axis=1, sort=False)

        self.coeffs = pd.concat([result_coeff_an, result_coeff_bn, result_coeff_cn, result_coeff_dn], axis=1,
                                sort=False)
        self.coeffs.index = wavelength_list

        self.cross_sections = pd.DataFrame([self.Qabs, self.Qext, self.Qscat]).T.rename(
            columns={0: 'Qabs', 1: 'Qext', 2: 'Qscat'})
        self.cross_sections.index = wavelength_list

    def compute_fields(self, wavelength, material, mp, points, components='EH', rat='Scatt'):
        # mp - multipoles

        self.Elec_flag = False
        self.Mag_flag = False
        self.full_field_flag = False

        if (components == 'EH'):
            self.Elec_flag = True
            self.Mag_flag = True
        elif (components == 'E'):
            self.Elec_flag = True
            self.Mag_flag = False
        elif (components == 'H'):
            self.Elec_flag = False
            self.Mag_flag = True
        else:
            raise ValueError('Wrong definition of the components.')

        if (rat == 'Scatt'):
            self.full_field_flag = False
        elif (rat == 'Full'):
            self.full_field_flag = True
        else:
            raise ValueError('Wrong definition of type of calculation. You must defined Scatt or Full.')

        self.check_parameters()

        #N_multipoles = self.N_multipoles
        m = self.m
        radius = self.a

        mp = np.arange(1, mp + 1, 1)
        mt = material.refractive_index(wavelength)

        k = 2 * np.pi * m / wavelength
        kt = 2 * np.pi * mt / wavelength

        alpha = k * radius
        beta = kt * radius

        an, bn, cn, dn = self.compute_coeffs(mt, m, mp, alpha, beta)

        def En(n, E0=1):
            return 1j ** n * (2 * n + 1) / (n * (n + 1)) * E0

        if (self.Elec_flag == True):
            Escatt_field = []

            for point in points:
                point_field = np.zeros(3)
                if (point[0] > radius):
                    for n in range(1, mp + 1):
                        point_field = point_field + (En(n) * (
                                    1j * np.multiply(an[n - 1], self.Ne1n(n, point, k)) - np.multiply(bn[n - 1],
                                                                                                      self.Mo1n(n,
                                                                                                                point,
                                                                                                                k))))
                    Escatt_field.append(point_field)

                else:
                    for n in range(1, mp + 1):
                        point_field = point_field + (En(n) * (
                                    np.multiply(cn[n - 1], self.Mo1n(n, point, kt, inside=True)) - 1j * np.multiply(
                                dn[n - 1], self.Ne1n(n, point, kt, inside=True))))
                    Escatt_field.append(point_field)

            Escatt_int_field = np.array(Escatt_field)

            self.E_scatt_int = Escatt_int_field

        if (self.Mag_flag == True):
            Hscatt_field = []

            for point in points:
                point_field = np.zeros(3)
                if (point[0] > radius):
                    for n in range(1, mp + 1):
                        point_field = point_field + (En(n) * (
                                    1j * np.multiply(bn[n - 1], self.No1n(n, point, k)) + np.multiply(an[n - 1], self.Me1n(n, point, k))))
                    Hscatt_field.append(k*self.c_const*point_field)

                else:
                    for n in range(1, mp + 1):
                        point_field = point_field + (En(n) * (
                                    np.multiply(dn[n - 1], self.Me1n(n, point, kt, inside=True)) + 1j * np.multiply(
                                cn[n - 1], self.No1n(n, point, kt, inside=True))))
                    Hscatt_field.append(-kt*self.c_const*point_field)

            Hscatt_int_field = np.array(Hscatt_field)
            self.H_scatt_int = Hscatt_int_field

        if (self.full_field_flag == True):
            E_pw_field = []

            for point in points:
                point_field = np.zeros(3)
                if (point[0] > radius):
                    for n in range(1, mp + 1):
                        point_field = point_field + (
                        (np.array(self.Mo1n(n, point, k)) - 1j * np.array(self.Ne1n(n, point, k))))
                    E_pw_field.append(point_field)
                else:
                    E_pw_field.append(np.array([0, 0, 0]))

            E_pw_field = np.array(E_pw_field)

            self.E_pw_field = E_pw_field

            H_pw_field = []

            for point in points:
                point_field = np.zeros(3)
                if (point[0] > radius):
                    for n in range(1, mp + 1):
                        point_field = point_field + (
                                    En(n) * (np.array(self.Me1n(n, point, k)) + 1j * np.array(self.No1n(n, point, k))))
                    H_pw_field.append(-k*self.c_const*point_field)
                else:
                    H_pw_field.append(np.array([0, 0, 0]))

            H_pw_field = np.array(H_pw_field)

            self.H_pw_field = H_pw_field

            self.E_total = self.E_pw_field + self.E_scatt_int
            self.H_total = self.H_pw_field + self.H_scatt_int

    def abs_fields(self):
        self.E_scatt_int_abs = np.einsum('ij,ij->i', self.E_scatt_int, np.conjugate(self.E_scatt_int))
        self.H_scatt_int_abs = np.einsum('ij,ij->i', self.H_scatt_int, np.conjugate(self.H_scatt_int))
        if (self.full_field_flag == True):
            self.E_total_abs = np.einsum('ij,ij->i', self.E_total, np.conjugate(self.E_total))
            self.H_total_abs = np.einsum('ij,ij->i', self.H_total, np.conjugate(self.H_total))

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

class MappingSpace(object):
    """
    How to use it:
    points = MappingSpace()
    points.mesh_gen()
    points.convert_cartesian_to_spherical()
    points.convert_spherical_to_cartesian()
    points.points

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    maps = MappingSpace()

    maps.mesh_gen(xmin = -1000, xmax= 1000, xstep=50, ymin = -1000, ymax = 1000, ystep = 50, plane = 'xy')
    ax.scatter(maps.points[:,0], maps.points[:,1], maps.points[:,2], marker = '.')
    maps.mesh_gen(xmin = -1000, xmax= 1000, xstep=50, ymin = -1000, ymax = 1000, ystep = 50, plane = 'xz')
    ax.scatter(maps.points[:,0], maps.points[:,1], maps.points[:,2], marker = '.')
    maps.mesh_gen(xmin = -1000, xmax= 1000, xstep=50, ymin = -1000, ymax = 1000, ystep = 50, plane = 'yz')
    ax.scatter(maps.points[:,0], maps.points[:,1], maps.points[:,2], marker = '.')
    """

    __version__ = 0.1

    # @staticmethod
    def mesh_gen(self, xmin=-500, xmax=500, xstep=10, ymin=-500, ymax=500, ystep=10, plane='xy'):
        """
        This function generates a rectangular plane with points in a regular mesh.
        parameters:
        return a 3xN vector with the points in the mesh.
        """

        xvalues = np.arange(xmin, xmax + xstep, xstep);
        yvalues = np.arange(ymin, ymax + ystep, ystep);
        xx, yy = np.meshgrid(xvalues, yvalues)

        if (plane == 'xy'):
            points_plane = np.vstack((xx.reshape(-1), yy.reshape(-1), 0 * yy.reshape(-1))).T
        elif (plane == 'xz'):
            points_plane = np.vstack((xx.reshape(-1), 0 * yy.reshape(-1), yy.reshape(-1))).T
        elif (plane == 'yz'):
            points_plane = np.vstack((0 * xx.reshape(-1), xx.reshape(-1), yy.reshape(-1))).T
        else:
            raise ValueError("The plane parameter is wrong. You must choose xy, xz or yz.")

        self.points = points_plane
        self.type = 'cartesian'

    def set_points(self, points, datatype):
        self.points = points
        if ((datatype == 'cartesian') or (datatype == 'spherical')):
            self.type = datatype
        else:
            raise ValueError('Datetype error. It must be cartesian or spherical.')

    # @staticmethod
    def convert_cartesian_to_spherical(self):
        if (self.type == 'cartesian'):
            # It is going to return (r, theta, phi) ISO convention
            ptsnew = np.hstack((self.points, np.zeros(self.points.shape)))
            xy = self.points[:, 0] ** 2 + self.points[:, 1] ** 2
            ptsnew[:, 3] = np.sqrt(xy + self.points[:, 2] ** 2)
            ptsnew[:, 4] = np.arctan2(np.sqrt(xy), self.points[:, 2])  # for elevation angle defined from Z-axis down
            # ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
            ptsnew[:, 5] = np.arctan2(self.points[:, 1], self.points[:, 0])
            self.points = ptsnew[:, 3:]
            self.type = 'spherical'
        else:
            print('The points are in spherical coordinates.')

    def convert_spherical_to_cartesian(self):
        if (self.type == 'spherical'):
            # Vector should be in the format (r, theta, phi) ISO convention
            X = self.points[:, 0] * np.sin(self.points[:, 1]) * np.cos(self.points[:, 2])
            Y = self.points[:, 0] * np.sin(self.points[:, 1]) * np.sin(self.points[:, 2])
            Z = self.points[:, 0] * np.cos(self.points[:, 1])
            self.points = np.vstack((X, Y, Z)).T
            self.type = 'cartesian'
        else:
            print('The points are in cartesian coordinates.')