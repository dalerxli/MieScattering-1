######################################################################
#
#
#
# Materials Class
#
# Mole Group
# December 2019
# Nuno de Sousa
######################################################################

import pandas as pd
import numpy as np

class Materials(object):
    """
    This is a general class that contains several functions that allow us to define and load the optical properties
    of materials.
    """

    material = None
    data = None
    path = None
    constant = False
    value = np.nan
    source = None

    def __init__(self, Symbol, path=None, value = None, source='refractiveindex.info'):
        """
        :param Symbol: Name that charaterizes the material
        :param path: path where the data is. Should be the address when using the refractiveindex or the path
        when we load a datafile. Data must be a csv with wavelength (column name wl), n (column name n) and k (column
        name k)
        :param value: It must contain a number when the source is "constant"
        :param source: It can be refractiveindex.info, datafile or constant.
        """
        self.material = Symbol
        self.path = path
        self.source = source
        if ((self.material != 'Si_palo') & (source == 'refractiveindex.info')):
            print('Downloading data from:')
            print(path)
            self.data = self.load_data_refractiveindex(path)
            print('Data Loaded.')

        if((self.material != 'Si_palo') & (source == 'datafile')):
            print('Downloading data from:')
            self.data = pd.read_csv(path).rename(columns = {'wl': 'wavelength'}).set_index('wavelength')
            print('Data Loaded.')

        if(self.material == 'Si_palo'):
            print('You have loaded the magic \'Si de palo\', i.e. wood silicon.')

        if((self.material != 'Si_palo') & (source == 'constant')):
            self.constant = True
            self.value = value


    def refractive_index(self, wavelength):
        """
        It returns the refractive index for a specific wavelength.
        :param wavelength: wavelength
        :return: refractive index (n + 1j*k)
        """
        if (self.material == 'Si_palo'):
            if ((wavelength >= 1E-6) & (wavelength <= 2E-6)):
                return 3.5 + 0.0j
            else:
                raise ValueError(
                    'It is not possible to compute the refractive index because the wavelength it is out of range.')

        elif(self.source == 'constant'):
            return self.value

        else:
            if ((wavelength >= self.data.index[0]) & (wavelength <= self.data.index[-1])):
                return self.interpolate(wavelength)
            else:
                raise ValueError(
                    'It is not possible to compute the refractive index because the wavelength it is out of range.')

    def load_data_refractiveindex(self, path):
        """
        :param path:
        :return:
        """
        data = pd.read_csv(path, header=None)
        n_index = data.index[data[1] == 'n'].values
        k_index = data.index[data[1] == 'k'].values

        # if database has n and k
        if ((np.size(n_index) != 0) & (np.size(k_index) != 0)):
            data_n = data[n_index[0]:k_index[0]]
            data_n.columns = data_n.iloc[0]
            data_n = data_n.drop(data_n.index[0])
            data_k = data[k_index[0]:]
            data_k.columns = data_k.iloc[0]
            data_k = data_k.drop(data_k.index[0])

            data_n = pd.DataFrame(data_n.reset_index()[['wl', 'n']].values, columns=['wavelength', 'n']).set_index(
                'wavelength')
            data_k = pd.DataFrame(data_k.reset_index()[['wl', 'k']].values, columns=['wavelength', 'k']).set_index(
                'wavelength')
            data_n['n'] = data_n['n'].astype(float)
            data_n.index = data_n.index.astype(float)
            data_k['k'] = data_k['k'].astype(float)
            data_k.index = data_k.index.astype(float)
            return pd.concat([data_n, data_k], axis=1, sort=False)

        # if database has only n
        if ((np.size(n_index) != 0) & (np.size(k_index) == 0)):
            data_n = data[n_index[0]:]
            data_n.columns = data_n.iloc[0]
            data_n = data_n.drop(data_n.index[0])
            data_n['n'] = data_n['n'].astype(float)
            data_n['wl'] = data_n['wl'].astype(float)
            data_n['k'] = 0
            data_n = data_n.rename(columns={'wl': 'wavelength'}).set_index('wavelength')
            return data_n

        # if database has only k
        if ((np.size(n_index) == 0) & (np.size(k_index) != 0)):
            data_k = data[k_index[0]:]
            data_k.columns = data_k.iloc[0]
            data_k = data_k.drop(data_k.index[0])
            data_k['k'] = data_k['k'].astype(float)
            data_k['wl'] = data_k['wl'].astype(float)
            data_k['n'] = 0
            data_k = data_k.rename(columns={'wl': 'wavelength'}).set_index('wavelength')
            return data_k[['n', 'k']]

    def interpolate(self, wavelength):
        """
        Linear interpolates the data to return the refractive index
        :param wavelength: wavelength to interpolate
        :return: refractive index (n + 1j*k)
        """
        return np.interp(wavelength, self.data.index.values, self.data['n'].values) + 1j * np.interp(wavelength,
                                                                                                     self.data.index.values,
                                                                                                     self.data[
                                                                                                         'k'].values)

    @staticmethod
    def polarization(volume=1, k=None, epsilon_part=1, epsilon_medium=1):

        # Clausius-Mossotti polarizability with radiative correction
        if (epsilon_part != -2):
            alpha0 = 3 * volume * (epsilon_part - 1 * epsilon_medium) / (epsilon_part + 2 * epsilon_medium)
            alpha = alpha0 / (1 - 1j * k ** 3 / (6 * np.pi) * alpha0)
            alpha0inv = 1 / (3 * volume * (epsilon_part - 1 * epsilon_medium) / (epsilon_part + 2 * epsilon_medium))

        # Resonant Particles
        elif (epsilon_part == -2):
            alpha = 1j * 6 * np.pi / k ** 3
            alpha0inv = 0

        else:
            raise ValueError('Something is wrong with the polarizability.')

        return alpha, alpha0inv