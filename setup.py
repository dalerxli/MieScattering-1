from setuptools import setup

setup(name='mole_mie',
      version='0.1.2',
      description='Scattering of a Spherical Particle using Mie Theory.',
      url='https://github.com/nunodsousa/MieScattering',
      author='Nuno de Sousa',
      author_email='nunodsousa@dipc.org',
      license='MIT',
      packages=['mole_mie'],
      install_requires=['termcolor', 'numpy', 'pandas', 'scipy>=1.3.1'],
      zip_safe=False)