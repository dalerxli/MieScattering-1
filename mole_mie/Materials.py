class Materials(object):
    material = None

    def __init__(self, Symbol):
        self.material = Symbol

    def refractive_index(self, wavelength):
        if (self.material == 'Si'):
            if ((wavelength >= 1E-6) & (wavelength <= 2E-6)):
                return 3.5 + 0.0j
            else:
                raise ValueError(
                    'It is not possible to compute the refractive index because the wavelength it is out of range.')
