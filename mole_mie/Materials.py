import pandas as pd
import numpy as np


class Materials(object):
    material = None
    data = None
    path = None

    def __init__(self, Symbol, path=None, source='refractiveindex.info'):
        self.material = Symbol
        self.path = path
        if ((self.material != 'Si_palo') & (source == 'refractiveindex.info')):
            print('Downloading data from:')
            print(path)
            self.data = self.load_data_refractiveindex(path)

    def refractive_index(self, wavelength):
        if (self.material == 'Si_palo'):
            if ((wavelength >= 1E-6) & (wavelength <= 2E-6)):
                return 3.5 + 0.0j
            else:
                raise ValueError(
                    'It is not possible to compute the refractive index because the wavelength it is out of range.')
        else:
            if ((wavelength >= self.data.index[0]) & (wavelength <= self.data.index[-1])):
                return self.interpolate(wavelength)
            else:
                raise ValueError(
                    'It is not possible to compute the refractive index because the wavelength it is out of range.')

    def load_data_refractiveindex(self, path):
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
        return np.interp(wavelength, self.data.index.values, self.data['n'].values) + 1j * np.interp(wavelength,
                                                                                                     self.data.index.values,
                                                                                                     self.data[
                                                                                                         'k'].values)
