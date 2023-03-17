from distutils.core import setup
from setuptools import find_packages, setup

setup(name='PiSlice',
      version='0.1',
      description='Estimate Pi and other population genomics statistics',
      url='https://github.com/ThomasBrazier/PiSlice',
      author='Thomas Brazier',
      author_email='thomas.brazier@univ-rennes1.fr',
      license='MIT',
      packages=['PiSlice'],
      extras_require=dict(tests=['pytest']),
      zip_safe=False)
